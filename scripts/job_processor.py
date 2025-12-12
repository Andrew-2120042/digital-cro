"""
Job Processor - Background worker for processing queued jobs

Run this script to process jobs from the queue:
    python scripts/job_processor.py

Or run in daemon mode:
    python scripts/job_processor.py --daemon
"""

import argparse
import time
import sys
from pathlib import Path
import logging

sys.path.insert(0, str(Path(__file__).parent.parent))

from scripts.utils.job_queue import JobQueue, JobStatus, JobType
from scripts.complete_workflow import run_complete_workflow
from scripts.utils.multi_target import run_multi_target_screening

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def process_single_target_job(job: dict, queue: JobQueue):
    """Process single-target screening job"""
    job_id = job['job_id']

    try:
        logger.info(f"Processing single-target job: {job_id}")

        # Update status to running
        queue.update_job_status(job_id, JobStatus.RUNNING.value)

        # Run workflow
        pdf_path = run_complete_workflow(
            pdb_id=job['target_proteins'][0],
            ligand_library_path=job['ligand_library_path'],
            output_dir=job['output_dir'],
            project_name=job['project_name'],
            client_name=job['client_name'],
            affinity_threshold=job['affinity_threshold'],
            max_workers=job['max_workers']
        )

        # Load results
        import pandas as pd
        results_csv = Path(job['output_dir']) / "final_results.csv"

        if results_csv.exists():
            df = pd.read_csv(results_csv)

            summary = {
                'total_molecules': len(df),
                'num_hits': len(df[df['binding_affinity'] <= job['affinity_threshold']]),
                'best_affinity': float(df['binding_affinity'].min()),
                'pdf_report': pdf_path
            }

            queue.update_job_status(
                job_id,
                JobStatus.COMPLETE.value,
                result_summary=summary
            )

            logger.info(f"✓ Job {job_id} complete: {summary['num_hits']} hits")
        else:
            raise Exception("Results file not found")

    except Exception as e:
        logger.error(f"✗ Job {job_id} failed: {e}")
        queue.update_job_status(
            job_id,
            JobStatus.FAILED.value,
            error_message=str(e)
        )


def process_multi_target_job(job: dict, queue: JobQueue):
    """Process multi-target screening job"""
    job_id = job['job_id']

    try:
        logger.info(f"Processing multi-target job: {job_id}")

        # Update status to running
        queue.update_job_status(job_id, JobStatus.RUNNING.value)

        # Run multi-target workflow
        results = run_multi_target_screening(
            target_proteins=job['target_proteins'],
            ligand_library_path=job['ligand_library_path'],
            output_dir=job['output_dir'],
            project_name=job['project_name'],
            client_name=job['client_name'],
            affinity_threshold=job['affinity_threshold'],
            max_workers_per_target=job['max_workers'],
            run_targets_parallel=False
        )

        # Extract summary
        summary = results['summary']

        queue.update_job_status(
            job_id,
            JobStatus.COMPLETE.value,
            result_summary=summary
        )

        logger.info(f"✓ Job {job_id} complete: {summary['successful_targets']} targets")

    except Exception as e:
        logger.error(f"✗ Job {job_id} failed: {e}")
        queue.update_job_status(
            job_id,
            JobStatus.FAILED.value,
            error_message=str(e)
        )


def process_next_job(queue: JobQueue) -> bool:
    """
    Process next job in queue.

    Returns:
        True if job was processed, False if queue empty
    """
    job = queue.get_next_queued_job()

    if not job:
        return False

    logger.info(f"Processing job: {job['job_id']} ({job['project_name']})")

    if job['job_type'] == JobType.SINGLE_TARGET.value:
        process_single_target_job(job, queue)
    elif job['job_type'] == JobType.MULTI_TARGET.value:
        process_multi_target_job(job, queue)
    else:
        logger.error(f"Unknown job type: {job['job_type']}")
        queue.update_job_status(
            job['job_id'],
            JobStatus.FAILED.value,
            error_message=f"Unknown job type: {job['job_type']}"
        )

    return True


def run_processor(daemon: bool = False, check_interval: int = 30):
    """
    Run job processor.

    Args:
        daemon: If True, run continuously
        check_interval: Seconds between queue checks (daemon mode)
    """
    queue = JobQueue()

    logger.info("Job processor started")
    logger.info(f"Daemon mode: {daemon}")

    if daemon:
        logger.info(f"Checking queue every {check_interval} seconds")

        while True:
            try:
                # Get queue stats
                stats = queue.get_queue_stats()
                logger.info(f"Queue: {stats['queued']} queued, {stats['running']} running, {stats['complete']} complete")

                # Process jobs while available
                while process_next_job(queue):
                    logger.info("Processed job, checking for more...")

                # Wait before next check
                time.sleep(check_interval)

            except KeyboardInterrupt:
                logger.info("Processor stopped by user")
                break
            except Exception as e:
                logger.error(f"Error in processor loop: {e}")
                time.sleep(check_interval)
    else:
        # Single-run mode: process all queued jobs
        processed = 0

        while process_next_job(queue):
            processed += 1

        logger.info(f"Processed {processed} jobs")

        # Show final stats
        stats = queue.get_queue_stats()
        logger.info(f"Final stats: {stats}")


def main():
    parser = argparse.ArgumentParser(
        description='Job Queue Processor for Digital CRO Platform'
    )

    parser.add_argument(
        '--daemon',
        action='store_true',
        help='Run in daemon mode (continuous processing)'
    )

    parser.add_argument(
        '--interval',
        type=int,
        default=30,
        help='Check interval in seconds (daemon mode)'
    )

    args = parser.parse_args()

    run_processor(daemon=args.daemon, check_interval=args.interval)


if __name__ == '__main__':
    main()
