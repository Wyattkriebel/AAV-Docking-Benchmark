#!/usr/bin/env python3
"""Audit archived inputs, freeze/prepare a separate8 A repair, or explicitly run it."""
from pathlib import Path
import argparse
import json
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from aav_benchmark.neighborhood_coverage import (
    AUDIT, audit_legacy, freeze_amendment, prepare_corrections, verify_corrections, run_corrections,
)
from aav_benchmark.provenance import write_json


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1])
    commands=parser.add_subparsers(dest='command',required=True)
    audit=commands.add_parser('audit');audit.add_argument('--check',action='store_true')
    commands.add_parser('freeze')
    commands.add_parser('prepare')
    commands.add_parser('verify')
    run=commands.add_parser('run');run.add_argument('--vina',type=Path,required=True);run.add_argument('--resume',action='store_true')
    args=parser.parse_args(argv)
    if args.command=='audit':
        result=audit_legacy(args.root);path=args.root/AUDIT
        if args.check:
            if result!=json.loads(path.read_text()):raise ValueError('Coverage audit is stale')
        else:
            if path.exists():raise FileExistsError('Use audit --check to preserve frozen source audit')
            write_json(path,result)
        print(json.dumps({k:result[k] for k in ['prepared_box_count','AAV_prepared_boxes','external_prepared_boxes','legacy_boxes_covering_8A','affected_box_keys']},indent=2))
    elif args.command=='freeze':
        result=freeze_amendment(args.root);print(json.dumps({'corrected_boxes':result['corrected_box_count'],'planned_runs':result['planned_run_count']}))
    elif args.command=='prepare':
        result=prepare_corrections(args.root);print(json.dumps({'status':result['status'],'corrected_boxes':result['corrected_box_count'],'expected_runs':result['expected_runs']}))
    elif args.command=='verify':print(json.dumps(verify_corrections(args.root),indent=2))
    else:
        result=run_corrections(args.root,args.vina,resume=args.resume);print(json.dumps({'boxes_executed':len(result)}))
    return 0


if __name__=='__main__':
    try:raise SystemExit(main())
    except (ValueError,OSError,KeyError,RuntimeError) as error:raise SystemExit(f'Neighborhood coverage failed: {error}')
