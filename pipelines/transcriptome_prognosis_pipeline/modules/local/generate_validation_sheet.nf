process GEN_VALIDATION_SHEET {
    tag 'gen_validation_sheet'

    publishDir "${params.run_dir}", mode: 'copy', overwrite: true

    input:
    val validation_ids_str
    val train_id

    output:
    path 'validation_sheet.csv', emit: validation_sheet

    script:
    """
     ${projectDir}/scripts/generate_validation_sheet.py \
      --rawdata-dir ${params.rawdata_dir} \
      --validation-ids '${validation_ids_str}' \
      --out validation_sheet.csv
    """
}
