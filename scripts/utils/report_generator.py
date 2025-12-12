"""
PDF Report Generator for Digital CRO.

Creates professional, client-ready PDF reports with:
    - Cover page with branding
    - Executive summary
    - Docking results with visualizations
    - ADMET analysis
    - Hit rankings
    - Methodology appendix
"""

from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER, TA_LEFT
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak,
    Table, TableStyle, KeepTogether
)
from reportlab.lib import colors
from datetime import datetime
from pathlib import Path
import pandas as pd
from typing import Optional, Dict, List
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class DrugDiscoveryReport:
    """
    Professional PDF report generator for drug discovery results.
    """

    def __init__(
        self,
        output_path: str,
        project_name: str = "Drug Discovery Project",
        client_name: str = "Client",
        page_size=letter
    ):
        """
        Initialize report generator.

        Args:
            output_path: Where to save PDF
            project_name: Name of the project
            client_name: Client/company name
            page_size: Page size (letter or A4)
        """
        self.output_path = Path(output_path)
        self.output_path.parent.mkdir(parents=True, exist_ok=True)

        self.project_name = project_name
        self.client_name = client_name
        self.page_size = page_size

        # Create PDF document
        self.doc = SimpleDocTemplate(
            str(self.output_path),
            pagesize=self.page_size,
            rightMargin=0.75*inch,
            leftMargin=0.75*inch,
            topMargin=0.75*inch,
            bottomMargin=0.75*inch
        )

        # Story (content to be added to PDF)
        self.story = []

        # Styles
        self.styles = getSampleStyleSheet()
        self._create_custom_styles()

        logger.info(f"Initialized report: {output_path}")

    def _create_custom_styles(self):
        """Create custom paragraph styles"""

        # Title style
        self.styles.add(ParagraphStyle(
            name='CustomTitle',
            parent=self.styles['Title'],
            fontSize=24,
            textColor=colors.HexColor('#1f4788'),
            spaceAfter=30,
            alignment=TA_CENTER
        ))

        # Heading1 style
        self.styles.add(ParagraphStyle(
            name='CustomHeading1',
            parent=self.styles['Heading1'],
            fontSize=16,
            textColor=colors.HexColor('#1f4788'),
            spaceAfter=12,
            spaceBefore=12
        ))

        # Heading2 style
        self.styles.add(ParagraphStyle(
            name='CustomHeading2',
            parent=self.styles['Heading2'],
            fontSize=14,
            textColor=colors.HexColor('#2e5c8a'),
            spaceAfter=10,
            spaceBefore=10
        ))

        # Body text
        self.styles.add(ParagraphStyle(
            name='Justify',
            parent=self.styles['BodyText'],
            alignment=TA_JUSTIFY,
            fontSize=11,
            leading=14
        ))

    def add_cover_page(
        self,
        target_protein: str,
        target_pdb_id: str,
        num_molecules_screened: int,
        date: str = None
    ):
        """
        Add cover page to report.

        Args:
            target_protein: Name of target protein
            target_pdb_id: PDB ID of target
            num_molecules_screened: Number of molecules screened
            date: Report date (default: today)
        """
        if date is None:
            date = datetime.now().strftime("%B %d, %Y")

        # Add logo/header space
        self.story.append(Spacer(1, 1.5*inch))

        # Title
        title = Paragraph(
            f"<b>Drug Discovery Report</b>",
            self.styles['CustomTitle']
        )
        self.story.append(title)
        self.story.append(Spacer(1, 0.3*inch))

        # Project name
        project = Paragraph(
            f"<b>{self.project_name}</b>",
            self.styles['CustomHeading1']
        )
        self.story.append(project)
        self.story.append(Spacer(1, 0.5*inch))

        # Project details
        details = [
            f"<b>Target Protein:</b> {target_protein} ({target_pdb_id})",
            f"<b>Molecules Screened:</b> {num_molecules_screened:,}",
            f"<b>Report Date:</b> {date}",
            f"<b>Prepared for:</b> {self.client_name}"
        ]

        for detail in details:
            p = Paragraph(detail, self.styles['Normal'])
            self.story.append(p)
            self.story.append(Spacer(1, 0.15*inch))

        self.story.append(PageBreak())

        logger.info("Added cover page")

    def add_executive_summary(
        self,
        summary_text: str,
        key_findings: List[str],
        top_hits: pd.DataFrame = None,
        num_top_hits: int = 5
    ):
        """
        Add executive summary section.

        Args:
            summary_text: Main summary paragraph
            key_findings: List of key findings (bullet points)
            top_hits: DataFrame with top hits (optional)
            num_top_hits: Number of top hits to highlight
        """
        # Section heading
        heading = Paragraph(
            "<b>Executive Summary</b>",
            self.styles['CustomHeading1']
        )
        self.story.append(heading)
        self.story.append(Spacer(1, 0.2*inch))

        # Summary text
        summary = Paragraph(summary_text, self.styles['Justify'])
        self.story.append(summary)
        self.story.append(Spacer(1, 0.2*inch))

        # Key findings
        findings_heading = Paragraph(
            "<b>Key Findings:</b>",
            self.styles['CustomHeading2']
        )
        self.story.append(findings_heading)
        self.story.append(Spacer(1, 0.1*inch))

        for finding in key_findings:
            bullet = Paragraph(f"• {finding}", self.styles['Normal'])
            self.story.append(bullet)
            self.story.append(Spacer(1, 0.05*inch))

        # Top hits table
        if top_hits is not None and not top_hits.empty:
            self.story.append(Spacer(1, 0.2*inch))

            top_heading = Paragraph(
                f"<b>Top {num_top_hits} Candidate Molecules:</b>",
                self.styles['CustomHeading2']
            )
            self.story.append(top_heading)
            self.story.append(Spacer(1, 0.1*inch))

            # Create table data
            table_data = [['Rank', 'Molecule', 'Binding Affinity', 'Drug-likeness', 'ADMET']]

            for idx, row in top_hits.head(num_top_hits).iterrows():
                rank = idx + 1
                mol_id = row.get('ligand_id', f'Molecule {rank}')
                affinity = f"{row.get('binding_affinity', 0):.2f} kcal/mol"
                qed = f"{row.get('qed_score', 0):.2f}" if 'qed_score' in row else 'N/A'

                # ADMET pass/fail
                lipinski_pass = row.get('lipinski_compliant', False)
                admet_status = '✓ Pass' if lipinski_pass else '✗ Fail'

                table_data.append([str(rank), mol_id, affinity, qed, admet_status])

            # Create table
            table = Table(table_data, colWidths=[0.6*inch, 1.5*inch, 1.3*inch, 1*inch, 1*inch])
            table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#1f4788')),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 10),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                ('GRID', (0, 0), (-1, -1), 1, colors.black),
                ('FONTSIZE', (0, 1), (-1, -1), 9),
            ]))

            self.story.append(table)

        self.story.append(PageBreak())

        logger.info("Added executive summary")

    def add_docking_results(
        self,
        df_results: pd.DataFrame,
        visualization_path: str = None,
        affinity_threshold: float = -7.0
    ):
        """
        Add docking results section.

        Args:
            df_results: DataFrame with docking results
            visualization_path: Path to molecule grid image
            affinity_threshold: Threshold for hit classification
        """
        # Section heading
        heading = Paragraph(
            "<b>Molecular Docking Results</b>",
            self.styles['CustomHeading1']
        )
        self.story.append(heading)
        self.story.append(Spacer(1, 0.2*inch))

        # Summary statistics
        total_molecules = len(df_results)
        hits = df_results[df_results['binding_affinity'] <= affinity_threshold]
        num_hits = len(hits)
        hit_rate = (num_hits / total_molecules) * 100 if total_molecules > 0 else 0

        best_affinity = df_results['binding_affinity'].min()
        mean_affinity = df_results['binding_affinity'].mean()

        summary_text = f"""
        A total of {total_molecules:,} molecules were screened against the target protein using
        molecular docking simulations. Binding affinities ranged from {best_affinity:.2f} to
        {df_results['binding_affinity'].max():.2f} kcal/mol. Using a threshold of {affinity_threshold} kcal/mol,
        {num_hits} molecules ({hit_rate:.1f}%) were classified as hits. The mean binding affinity
        was {mean_affinity:.2f} kcal/mol.
        """

        summary = Paragraph(summary_text, self.styles['Justify'])
        self.story.append(summary)
        self.story.append(Spacer(1, 0.2*inch))

        # Add visualization if provided
        if visualization_path and Path(visualization_path).exists():
            img = Image(visualization_path, width=6*inch, height=4*inch)
            self.story.append(img)
            self.story.append(Spacer(1, 0.1*inch))

            caption = Paragraph(
                "<i>Figure 1: Top-ranked molecules with binding affinities</i>",
                self.styles['Normal']
            )
            self.story.append(caption)

        self.story.append(PageBreak())

        logger.info("Added docking results section")

    def add_admet_analysis(
        self,
        df_results: pd.DataFrame,
        admet_summary_path: str = None,
        radar_plot_path: str = None
    ):
        """
        Add ADMET analysis section.

        Args:
            df_results: DataFrame with ADMET predictions
            admet_summary_path: Path to ADMET summary visualization
            radar_plot_path: Path to radar plot for top molecule
        """
        # Section heading
        heading = Paragraph(
            "<b>ADMET Analysis</b>",
            self.styles['CustomHeading1']
        )
        self.story.append(heading)
        self.story.append(Spacer(1, 0.2*inch))

        # ADMET summary
        if 'lipinski_compliant' in df_results.columns:
            lipinski_pass = df_results['lipinski_compliant'].sum()
            lipinski_rate = (lipinski_pass / len(df_results)) * 100

            mean_qed = df_results.get('qed_score', pd.Series([0])).mean()

            bbb_count = df_results.get('bbb_penetrant', pd.Series([False])).sum()
            bbb_rate = (bbb_count / len(df_results)) * 100

            high_bioavail = (df_results.get('oral_bioavailability', pd.Series([''])) == 'High').sum()
            bioavail_rate = (high_bioavail / len(df_results)) * 100

            summary_text = f"""
            Drug-likeness assessment was performed on all screened molecules using Lipinski's Rule of Five
            and quantitative estimate of drug-likeness (QED). {lipinski_pass} molecules ({lipinski_rate:.1f}%)
            passed Lipinski's criteria, with a mean QED score of {mean_qed:.2f}. Blood-brain barrier (BBB)
            penetration was predicted for {bbb_count} molecules ({bbb_rate:.1f}%), and {high_bioavail}
            molecules ({bioavail_rate:.1f}%) were predicted to have high oral bioavailability.
            """

            summary = Paragraph(summary_text, self.styles['Justify'])
            self.story.append(summary)
            self.story.append(Spacer(1, 0.2*inch))

        # Add ADMET summary visualization
        if admet_summary_path and Path(admet_summary_path).exists():
            img = Image(admet_summary_path, width=6.5*inch, height=4.5*inch)
            self.story.append(img)
            self.story.append(Spacer(1, 0.1*inch))

            caption = Paragraph(
                "<i>Figure 2: ADMET property distributions and compliance rates</i>",
                self.styles['Normal']
            )
            self.story.append(caption)
            self.story.append(Spacer(1, 0.2*inch))

        # Add radar plot for top molecule
        if radar_plot_path and Path(radar_plot_path).exists():
            img = Image(radar_plot_path, width=4*inch, height=4*inch)
            self.story.append(img)
            self.story.append(Spacer(1, 0.1*inch))

            caption = Paragraph(
                "<i>Figure 3: ADMET profile of top-ranked molecule</i>",
                self.styles['Normal']
            )
            self.story.append(caption)

        self.story.append(PageBreak())

        logger.info("Added ADMET analysis section")

    def add_methodology(
        self,
        docking_params: Dict = None
    ):
        """
        Add methodology appendix.

        Args:
            docking_params: Dictionary with docking parameters
        """
        # Section heading
        heading = Paragraph(
            "<b>Methodology</b>",
            self.styles['CustomHeading1']
        )
        self.story.append(heading)
        self.story.append(Spacer(1, 0.2*inch))

        # Molecular Docking
        subheading = Paragraph(
            "<b>Molecular Docking</b>",
            self.styles['CustomHeading2']
        )
        self.story.append(subheading)
        self.story.append(Spacer(1, 0.1*inch))

        docking_text = """
        Molecular docking simulations were performed using AutoDock Vina. The target protein structure
        was prepared by removing water molecules, adding hydrogen atoms, and assigning partial charges.
        The binding pocket was identified using ligand-based pocket detection. All ligands were converted
        to PDBQT format with energy minimization and partial charge assignment using Gasteiger method.
        Docking was performed with exhaustiveness of 8, generating 9 binding modes per ligand.
        """

        text = Paragraph(docking_text, self.styles['Justify'])
        self.story.append(text)
        self.story.append(Spacer(1, 0.2*inch))

        # ADMET Predictions
        subheading = Paragraph(
            "<b>ADMET Predictions</b>",
            self.styles['CustomHeading2']
        )
        self.story.append(subheading)
        self.story.append(Spacer(1, 0.1*inch))

        admet_text = """
        Drug-likeness was assessed using Lipinski's Rule of Five (molecular weight ≤ 500 Da, LogP ≤ 5,
        hydrogen bond donors ≤ 5, hydrogen bond acceptors ≤ 10). Quantitative estimate of drug-likeness
        (QED) was calculated using RDKit's implementation. Blood-brain barrier penetration was predicted
        using a multi-criteria model (TPSA < 90 Ų, MW < 450 Da, optimal LogP 1-3). Oral bioavailability
        was estimated using combined Lipinski and Veber criteria (rotatable bonds ≤ 10, TPSA ≤ 140 Ų).
        Synthetic accessibility scores were calculated based on molecular complexity.
        """

        text = Paragraph(admet_text, self.styles['Justify'])
        self.story.append(text)

        logger.info("Added methodology section")

    def generate(self):
        """Build and save the PDF report"""
        try:
            self.doc.build(self.story)
            logger.info(f"✓ Report generated successfully: {self.output_path}")
            return str(self.output_path)
        except Exception as e:
            logger.error(f"Error generating report: {e}")
            raise


def generate_complete_report(
    df_results: pd.DataFrame,
    output_path: str,
    project_name: str,
    client_name: str,
    target_protein: str,
    target_pdb_id: str,
    visualization_dir: str = None,
    affinity_threshold: float = -7.0
) -> str:
    """
    Generate complete drug discovery report.

    Args:
        df_results: DataFrame with docking + ADMET results
        output_path: Where to save PDF
        project_name: Project name
        client_name: Client name
        target_protein: Target protein name
        target_pdb_id: PDB ID
        visualization_dir: Directory with visualization images
        affinity_threshold: Threshold for hit classification

    Returns:
        Path to generated PDF
    """
    # Create report
    report = DrugDiscoveryReport(
        output_path=output_path,
        project_name=project_name,
        client_name=client_name
    )

    # Cover page
    report.add_cover_page(
        target_protein=target_protein,
        target_pdb_id=target_pdb_id,
        num_molecules_screened=len(df_results)
    )

    # Executive summary
    hits = df_results[df_results['binding_affinity'] <= affinity_threshold]
    num_hits = len(hits)
    best_affinity = df_results['binding_affinity'].min()

    summary_text = f"""
    This report presents the results of a computational drug discovery campaign targeting {target_protein}.
    A virtual screening of {len(df_results):,} drug-like molecules identified {num_hits} promising candidates
    with strong predicted binding affinity (≤ {affinity_threshold} kcal/mol). The top-ranked molecule showed
    a binding affinity of {best_affinity:.2f} kcal/mol. All hits were further evaluated for drug-likeness,
    ADMET properties, and synthetic accessibility.
    """

    key_findings = [
        f"{num_hits} molecules passed binding affinity threshold ({affinity_threshold} kcal/mol)",
        f"Best binding affinity: {best_affinity:.2f} kcal/mol",
        f"Mean QED score: {df_results.get('qed_score', pd.Series([0])).mean():.2f}",
        f"{df_results.get('lipinski_compliant', pd.Series([False])).sum()} molecules passed Lipinski's Rule of Five"
    ]

    report.add_executive_summary(
        summary_text=summary_text,
        key_findings=key_findings,
        top_hits=df_results.sort_values('binding_affinity').head(10)
    )

    # Docking results
    viz_path = None
    if visualization_dir:
        viz_dir = Path(visualization_dir)
        potential_paths = [
            viz_dir / 'top_hits.png',
            viz_dir / 'molecules_grid.png',
            viz_dir / 'docking_results.png'
        ]
        for p in potential_paths:
            if p.exists():
                viz_path = str(p)
                break

    report.add_docking_results(
        df_results=df_results,
        visualization_path=viz_path,
        affinity_threshold=affinity_threshold
    )

    # ADMET analysis
    admet_summary_path = None
    radar_path = None
    if visualization_dir:
        viz_dir = Path(visualization_dir)
        if (viz_dir / 'admet_summary.png').exists():
            admet_summary_path = str(viz_dir / 'admet_summary.png')

        # Look for any radar plot
        radar_files = list(viz_dir.glob('admet_radar_*.png'))
        if radar_files:
            radar_path = str(radar_files[0])

    report.add_admet_analysis(
        df_results=df_results,
        admet_summary_path=admet_summary_path,
        radar_plot_path=radar_path
    )

    # Methodology
    report.add_methodology()

    # Generate PDF
    return report.generate()
