from __future__ import annotations
import os, re

# NCBI / Eutils
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
NCBI_API_KEY = os.getenv("NCBI_API_KEY", "").strip()
TOOL_NAME = os.getenv("NCBI_TOOL", "urbanscope-srr-harvester")
NCBI_EMAIL = os.getenv("NCBI_EMAIL", "")

"""
DEFAULT_QUERY = (
    '((urban OR city OR cities OR subway OR transit OR "built environment" '
    'OR wastewater OR sewage OR stormwater OR "public transit" OR "surface swab") '
    'AND (microbiome OR metagenom* OR metatranscriptom* OR "shotgun metagenomic" '
    'OR "environmental swab" OR "RNA-seq" OR "amplicon"))'
)
"""

"""
DEFAULT_QUERY = (
    '('
    '('
        'urban OR city OR cities OR metropolitan OR municipal OR '
        '"built environment" OR infrastructure OR '
        'subway OR metro OR transit OR railway OR airport OR '
        'wastewater OR sewage OR stormwater OR '
        'sidewalk OR street OR road OR pavement OR '
        'building OR housing OR hospital OR school OR '
        'surface OR fomite OR air OR aerosol'
    ') '
    'AND '
    '('
        'metagenom* OR metatranscriptom* OR '
        '"shotgun sequencing" OR '
        '"shotgun metagenomic" OR '
        '"environmental sequencing" OR '
        '"environmental DNA" OR eDNA OR '
        '"RNA sequencing" OR RNA-seq OR '
        'microbiome OR microbiota OR '
        'amplicon OR 16S OR 18S OR ITS'
    ')'
    ')'
)
"""

"""
DEFAULT_QUERY = (
    '('
      '('
        # Urban/city constraint (strong)
        '(urban OR "built environment" OR city OR cities OR metropolitan OR municipal OR '
        'subway OR metro OR transit OR railway OR airport OR '
        'wastewater OR sewage OR stormwater OR '
        '"public transport" OR "public transit" OR '
        'sidewalk OR street OR pavement OR '
        'building OR housing OR hospital OR school OR '
        '"surface swab" OR surface OR fomite OR air OR aerosol)'
      ') '
      'AND '
      '('
        # Shotgun WGS metagenomics OR metatranscriptomics constraint (strong)
        '("shotgun metagenom*" OR "shotgun sequencing" OR "whole genome shotgun" OR WGS OR '
        'metagenom* OR metatranscriptom* OR "total RNA" OR "metatranscriptome")'
      ') '
      'AND '
      '('
        # Environmental sampling signal (helps keep it non-host and non-clinical)
        '(environment* OR "built environment" OR wastewater OR sewage OR stormwater OR '
        'surface OR swab OR air OR aerosol)'
      ') '
      'NOT '
      '('
        # Exclude amplicon / marker-gene / targeted assays
        'amplicon OR "marker gene" OR "16S" OR "16S rRNA" OR "18S" OR "ITS" OR '
        '"COI" OR "barcod*" OR "qPCR" OR "PCR amplicon" OR "V3-V4" OR "V4 region"'
      ')'
    ')'
)

"""


DEFAULT_QUERY = (
    '('
      # -------------------------------
      # 1) MUST be urban / city context
      # -------------------------------
      '('
        '"urban" OR "city" OR "cities" OR metropolitan OR municipal OR '
        '"built environment" OR '
        'subway OR metro OR transit OR railway OR airport OR '
        '"public transit" OR "public transport" OR '
        'wastewater OR sewage OR stormwater OR '
        'street OR sidewalk OR pavement OR '
        'building OR buildings OR housing OR '
        '"surface swab" OR fomite OR '
        'air OR aerosol'
      ') '
      'AND '

      # ------------------------------------------------
      # 2) MUST be shotgun metagenomics / metatranscriptomics
      # ------------------------------------------------
      '('
        '"whole genome shotgun" OR '
        '"shotgun metagenom*" OR '
        '"shotgun sequencing" OR '
        'metagenom* OR '
        'metatranscriptom* OR '
        '"total RNA sequencing" OR '
        '"metatranscriptome sequencing"'
      ') '
      'AND '

      # ------------------------------------------------
      # 3) MUST be environmental (not host clinical RNA-seq)
      # ------------------------------------------------
      '('
        'environment* OR '
        '"built environment" OR '
        'wastewater OR sewage OR stormwater OR '
        'surface OR swab OR '
        'air OR aerosol'
      ') '

      # -------------------------------
      # 4) HARD EXCLUSIONS
      # -------------------------------
      'NOT '
      '('
        # Marker gene / amplicon
        'amplicon OR '
        '"marker gene" OR '
        '"16S" OR "16S rRNA" OR '
        '"18S" OR '
        '"ITS" OR '
        '"V3-V4" OR "V4 region" OR '
        '"barcod*" OR '

        # Host-focused RNA-seq / genomics
        '"RNA-seq" OR '
        '"single-cell" OR '
        'scRNA OR '
        '"whole exome" OR '
        'WES OR '

        # Common non-urban environments
        'soil OR '
        'sediment OR '
        'marine OR '
        'ocean OR '
        'freshwater OR '
        'river OR '
        'lake OR '
        'forest OR '
        'agricultur* OR '
        'farm OR '
        'plant OR '
        'rhizosphere'
      ')'
    ')'
)







# Paths
DATA_DIR = "data"
DOCS_DIR = "docs"
DB_DIR = f"{DOCS_DIR}/db"
CACHE_DIR = f"{DATA_DIR}/cache"
DEBUG_DIR = f"{DATA_DIR}/debug"
DOCS_DEBUG_DIR = f"{DOCS_DIR}/debug"

SEEN_SRA_UIDS = f"{DATA_DIR}/seen_sra_uids.txt"
SEEN_SRR_RUNS = f"{DATA_DIR}/seen_srr_runs.txt"

BIOSAMPLE_CACHE = f"{CACHE_DIR}/biosample.json"
BIOPROJECT_CACHE = f"{CACHE_DIR}/bioproject.json"
BIOPROJECT_UID_CACHE = f"{CACHE_DIR}/bioproject_uid.json"

DOCS_LATEST_SRR = f"{DOCS_DIR}/latest_srr.json"
DOCS_LATEST_DEBUG = f"{DOCS_DEBUG_DIR}/latest_report.json"

BIOPROJECT_RE = re.compile(r"\bPRJ(?:NA|EB|DB)\d+\b", re.I)

# 50MB limit
MAX_OUTPUT_BYTES = 50 * 1024 * 1024
