#!/usr/bin/env bash
set -euo pipefail

THREADS="${THREADS:-8}"
MAX_SIZE="${MAX_SIZE:-20G}"
RUN_LIST="${RUN_LIST:-config/sra_runs.txt}"

SRA_ROOT="${SRA_ROOT:-data/sra}"
SRA_CACHE="${SRA_ROOT}/cache"
FASTQ_DIR="${SRA_ROOT}/fastq"
LOG_DIR="${SRA_ROOT}/logs"

mkdir -p "${SRA_CACHE}" "${FASTQ_DIR}" "${LOG_DIR}"

while read -r RUN; do
  [[ -z "${RUN}" || "${RUN}" =~ ^# ]] && continue

  echo "Downloading ${RUN}"
  prefetch \
    --max-size "${MAX_SIZE}" \
    --output-directory "${SRA_CACHE}" \
    "${RUN}" \
    2>&1 | tee "${LOG_DIR}/${RUN}.prefetch.log"

  echo "Converting ${RUN} to paired FASTQ"
  fasterq-dump \
    "${SRA_CACHE}/${RUN}/${RUN}.sra" \
    --split-files \
    --threads "${THREADS}" \
    --outdir "${FASTQ_DIR}" \
    --temp /tmp \
    2>&1 | tee "${LOG_DIR}/${RUN}.fasterq-dump.log"

  echo "Compressing ${RUN} FASTQ files"
  pigz -p "${THREADS}" "${FASTQ_DIR}/${RUN}"_*.fastq
done < "${RUN_LIST}"

echo "Downloaded and converted all SRA runs in ${RUN_LIST}"
