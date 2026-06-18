process RAW_DATA_COPY {
    tag 'copy_rawdata'

    publishDir "${params.rawdata_dir}", mode: 'copy', overwrite: false

    input:
    val(clean_root)
    val(train_id)
    val(validation_ids_str)

    output:
    path '01.*_count.csv',          optional: true
    path '02.*_group.csv',          optional: true
    path '03.*_tpmlog2.csv',        optional: true
    path '04.*_fpkmlog2.csv',       optional: true
    path '05.*_survival.csv',       optional: true
    path '06.*_clinical.csv',       optional: true
    path '*/07.*_expr.csv',         optional: true
    path '*/09.*_survival.csv',     optional: true

    script:
    """
    # Copy training data
    train_src="${clean_root}/Train_set/${train_id}"
    if [ -d "\$train_src" ]; then
        cp "\$train_src"/*.csv . 2>/dev/null || true
    fi

    # Copy validation cohorts
    IFS=',' read -ra VIDS <<< "${validation_ids_str}"
    for vid in "\${VIDS[@]}"; do
        vid=\$(echo "\$vid" | xargs)
        [ -z "\$vid" ] && continue
        src="${clean_root}/Validation_set/\${vid}"
        if [ -d "\$src" ]; then
            mkdir -p "\$vid"
            cp "\$src"/*.csv "\$vid/" 2>/dev/null || true
        fi
    done
    """
}
