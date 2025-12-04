"""
Job Queue System for Digital CRO Platform

Simple FIFO queue for batch screening jobs.
Uses JSON file storage (can upgrade to SQLite/Postgres later).
"""

import json
import uuid
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Optional
import logging
from dataclasses import dataclass, asdict
from enum import Enum

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class JobStatus(Enum):
    QUEUED = "queued"
    RUNNING = "running"
    COMPLETE = "complete"
    FAILED = "failed"
    CANCELLED = "cancelled"


class JobType(Enum):
    SINGLE_TARGET = "single_target"
    MULTI_TARGET = "multi_target"


@dataclass
class Job:
    """Screening job"""
    job_id: str
    job_type: str
    status: str
    created_at: str
    started_at: Optional[str] = None
    completed_at: Optional[str] = None

    # Job parameters
    project_name: str = ""
    client_name: str = ""
    target_proteins: List[str] = None
    ligand_library_path: str = ""
    output_dir: str = ""
    affinity_threshold: float = -7.0
    max_workers: int = 4
    use_consensus: bool = False

    # Results
    result_summary: Optional[Dict] = None
    error_message: Optional[str] = None

    def __post_init__(self):
        if self.target_proteins is None:
            self.target_proteins = []


class JobQueue:
    """
    Simple job queue manager.

    Uses JSON file for persistence (scalable to 1000s of jobs).
    """

    def __init__(self, queue_dir: str = "data/job_queue"):
        self.queue_dir = Path(queue_dir)
        self.queue_dir.mkdir(parents=True, exist_ok=True)

        self.jobs_file = self.queue_dir / "jobs.json"
        self.lock_file = self.queue_dir / ".lock"

        # Initialize jobs file if doesn't exist
        if not self.jobs_file.exists():
            self._save_jobs({})

    def _load_jobs(self) -> Dict[str, Dict]:
        """Load all jobs from storage"""
        try:
            with open(self.jobs_file, 'r') as f:
                return json.load(f)
        except:
            return {}

    def _save_jobs(self, jobs: Dict[str, Dict]):
        """Save all jobs to storage"""
        with open(self.jobs_file, 'w') as f:
            json.dump(jobs, f, indent=2)

    def submit_single_target_job(
        self,
        pdb_id: str,
        ligand_library_path: str,
        project_name: str = "Screening Job",
        client_name: str = "Client",
        affinity_threshold: float = -7.0,
        max_workers: int = 4,
        use_consensus: bool = False
    ) -> str:
        """
        Submit single-target screening job.

        Returns:
            job_id: Unique identifier for this job
        """
        job_id = str(uuid.uuid4())

        output_dir = self.queue_dir.parent / "outputs" / f"job_{job_id}"

        job = Job(
            job_id=job_id,
            job_type=JobType.SINGLE_TARGET.value,
            status=JobStatus.QUEUED.value,
            created_at=datetime.now().isoformat(),
            project_name=project_name,
            client_name=client_name,
            target_proteins=[pdb_id],
            ligand_library_path=ligand_library_path,
            output_dir=str(output_dir),
            affinity_threshold=affinity_threshold,
            max_workers=max_workers,
            use_consensus=use_consensus
        )

        # Save to queue
        jobs = self._load_jobs()
        jobs[job_id] = asdict(job)
        self._save_jobs(jobs)

        logger.info(f"✓ Job submitted: {job_id} ({project_name})")

        return job_id

    def submit_multi_target_job(
        self,
        target_proteins: List[str],
        ligand_library_path: str,
        project_name: str = "Multi-Target Job",
        client_name: str = "Client",
        affinity_threshold: float = -7.0,
        max_workers: int = 4,
        use_consensus: bool = False
    ) -> str:
        """
        Submit multi-target screening job.

        Returns:
            job_id: Unique identifier for this job
        """
        job_id = str(uuid.uuid4())

        output_dir = self.queue_dir.parent / "outputs" / f"job_{job_id}"

        job = Job(
            job_id=job_id,
            job_type=JobType.MULTI_TARGET.value,
            status=JobStatus.QUEUED.value,
            created_at=datetime.now().isoformat(),
            project_name=project_name,
            client_name=client_name,
            target_proteins=target_proteins,
            ligand_library_path=ligand_library_path,
            output_dir=str(output_dir),
            affinity_threshold=affinity_threshold,
            max_workers=max_workers,
            use_consensus=use_consensus
        )

        # Save to queue
        jobs = self._load_jobs()
        jobs[job_id] = asdict(job)
        self._save_jobs(jobs)

        logger.info(f"✓ Multi-target job submitted: {job_id} ({project_name})")

        return job_id

    def get_job(self, job_id: str) -> Optional[Dict]:
        """Get job by ID"""
        jobs = self._load_jobs()
        return jobs.get(job_id)

    def list_jobs(
        self,
        status: Optional[str] = None,
        limit: int = 100
    ) -> List[Dict]:
        """
        List all jobs, optionally filtered by status.

        Returns jobs sorted by created_at (newest first).
        """
        jobs = self._load_jobs()
        job_list = list(jobs.values())

        # Filter by status
        if status:
            job_list = [j for j in job_list if j['status'] == status]

        # Sort by created_at (newest first)
        job_list.sort(key=lambda x: x['created_at'], reverse=True)

        return job_list[:limit]

    def get_next_queued_job(self) -> Optional[Dict]:
        """Get next job to process (FIFO)"""
        jobs = self._load_jobs()

        queued_jobs = [
            j for j in jobs.values()
            if j['status'] == JobStatus.QUEUED.value
        ]

        if not queued_jobs:
            return None

        # Sort by created_at (oldest first - FIFO)
        queued_jobs.sort(key=lambda x: x['created_at'])

        return queued_jobs[0]

    def update_job_status(
        self,
        job_id: str,
        status: str,
        result_summary: Optional[Dict] = None,
        error_message: Optional[str] = None
    ):
        """Update job status and results"""
        jobs = self._load_jobs()

        if job_id not in jobs:
            raise ValueError(f"Job {job_id} not found")

        job = jobs[job_id]
        job['status'] = status

        if status == JobStatus.RUNNING.value and not job.get('started_at'):
            job['started_at'] = datetime.now().isoformat()

        if status in [JobStatus.COMPLETE.value, JobStatus.FAILED.value]:
            job['completed_at'] = datetime.now().isoformat()

        if result_summary:
            job['result_summary'] = result_summary

        if error_message:
            job['error_message'] = error_message

        self._save_jobs(jobs)

        logger.info(f"✓ Job {job_id} status: {status}")

    def cancel_job(self, job_id: str):
        """Cancel a queued job"""
        job = self.get_job(job_id)

        if not job:
            raise ValueError(f"Job {job_id} not found")

        if job['status'] not in [JobStatus.QUEUED.value, JobStatus.RUNNING.value]:
            raise ValueError(f"Cannot cancel job in status: {job['status']}")

        self.update_job_status(job_id, JobStatus.CANCELLED.value)

        logger.info(f"✓ Job {job_id} cancelled")

    def delete_job(self, job_id: str):
        """Delete a job (and its results)"""
        jobs = self._load_jobs()

        if job_id not in jobs:
            raise ValueError(f"Job {job_id} not found")

        # Delete output directory
        job = jobs[job_id]
        output_dir = Path(job['output_dir'])

        if output_dir.exists():
            import shutil
            shutil.rmtree(output_dir)
            logger.info(f"✓ Deleted output: {output_dir}")

        # Remove from queue
        del jobs[job_id]
        self._save_jobs(jobs)

        logger.info(f"✓ Job {job_id} deleted")

    def get_queue_stats(self) -> Dict:
        """Get queue statistics"""
        jobs = self._load_jobs()

        stats = {
            'total_jobs': len(jobs),
            'queued': sum(1 for j in jobs.values() if j['status'] == JobStatus.QUEUED.value),
            'running': sum(1 for j in jobs.values() if j['status'] == JobStatus.RUNNING.value),
            'complete': sum(1 for j in jobs.values() if j['status'] == JobStatus.COMPLETE.value),
            'failed': sum(1 for j in jobs.values() if j['status'] == JobStatus.FAILED.value),
            'cancelled': sum(1 for j in jobs.values() if j['status'] == JobStatus.CANCELLED.value)
        }

        return stats
