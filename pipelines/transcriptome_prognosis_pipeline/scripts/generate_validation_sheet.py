#!/usr/bin/env python3
"""Auto-generate validation_sheet.csv based on survival time ranges."""
import csv, glob, os, sys, argparse

def read_surv_range(surv_file):
    """Read survival file, return (min_OS.time, max_OS.time)."""
    if not os.path.exists(surv_file):
        return None, None
    try:
        with open(surv_file) as f:
            reader = csv.DictReader(f)
            times = []
            for row in reader:
                t = row.get('OS.time', row.get('OS_TIME', ''))
                if t:
                    try: times.append(float(t))
                    except ValueError: pass
            if times:
                return min(times), max(times)
    except Exception:
        pass
    return None, None

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--rawdata-dir', required=True)
    p.add_argument('--validation-ids', required=True, help='comma-separated cohort IDs')
    p.add_argument('--out', required=True, help='output CSV path')
    args = p.parse_args()

    rawdata = args.rawdata_dir
    valid_ids = [v.strip() for v in args.validation_ids.split(',') if v.strip()]

    # Find training survival file
    train_surv = None
    for f in glob.glob(os.path.join(rawdata, '05.*_survival.csv')):
        train_surv = f; break
    if not train_surv:
        print("ERROR: training survival file not found")
        sys.exit(1)

    t_min, t_max = read_surv_range(train_surv)
    if t_min is None:
        print("ERROR: cannot parse training survival")
        sys.exit(1)

    rows = []
    for cid in valid_ids:
        cdir = os.path.join(rawdata, cid)
        if not os.path.isdir(cdir):
            print(f"WARN: validation dir not found: {cdir}, skipping")
            continue

        # Find cohort expression/survival files
        expr_f = surv_f = None
        for f in glob.glob(os.path.join(cdir, '07.*_expr.csv')):
            expr_f = f; break
        for f in glob.glob(os.path.join(cdir, '09.*_survival.csv')):
            surv_f = f; break
        if not expr_f or not surv_f:
            print(f"WARN: missing files for {cid}, skipping")
            continue

        v_min, v_max = read_surv_range(surv_f)
        if v_min is None:
            print(f"WARN: cannot parse survival for {cid}, skipping")
            continue

        # Rule: exclude if valid max < 5yr AND valid min > 1yr
        if v_max < 1825 and v_min > 365:
            print(f"EXCLUDE {cid}: max={v_max} < 5yr AND min={v_min} > 1yr")
            continue

        need_135 = False
        need_357 = False

        if t_max < 1825:
            # Training max < 5yr
            if v_min <= 365:
                need_135 = True
        elif t_min > 365:
            # Training min > 1yr
            if v_max >= 1825:
                need_357 = True
        else:
            # Training spans both
            if v_max < 1825:
                need_135 = True
            elif v_min > 365:
                need_357 = True
            else:
                need_135 = True
                need_357 = True

        if not need_135 and not need_357:
            print(f"SKIP {cid}: no applicable time set")
            continue

        if need_135:
            rows.append([cid, os.path.abspath(expr_f), os.path.abspath(surv_f), ''])
            print(f"ADD {cid} 135y")
        if need_357:
            rid = f"{cid}_357y" if not need_135 else cid
            td = '1095;1825;2555' if need_135 else '1095;1825;2555'
            if need_135 and need_357:
                rid = f"{cid}_357y"
            rows.append([rid, os.path.abspath(expr_f), os.path.abspath(surv_f),
                        '1095;1825;2555' if not need_135 else ('1095;1825;2555' if need_357 else '')])
            if need_135 and need_357:
                rows[-1][3] = '1095;1825;2555'
                print(f"ADD {cid} 357y as {rid}")
            elif need_357:
                rows[-1][3] = '1095;1825;2555'
                print(f"ADD {cid} 357y")

    with open(args.out, 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['cohort_id', 'expr_file', 'surv_file', 'time_roc_days'])
        for row in rows:
            w.writerow(row)
    print(f"Wrote {len(rows)} rows to {args.out}")

if __name__ == '__main__':
    main()
