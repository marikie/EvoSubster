"""Explicit identity conditions for PCA; JSON request/response bridge for R."""
import json
from pathlib import Path
import sys
from decimal import Decimal
from train_identity import IDENTITY_SOURCE
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'select'))
from trio_filter import extract_genus_label, identity_gaps, validate_identity_threshold


def evaluate(request):
    args = request['options']
    common = args.get('idt_threshold')
    thresholds = {
        key: args.get(key) if args.get(key) is not None else common
        for key in ('ingroup_idt_threshold', 'outgroup_idt_threshold')
    }
    thresholds['min_idt_gap_between_out_in'] = args.get('min_idt_gap_between_out_in')
    for value in [common, *thresholds.values()]:
        if value is not None:
            validate_identity_threshold(value)
    values = {}
    for key, value in (request.get('identities') or {}).items():
        try:
            values[key] = validate_identity_threshold(value)
        except (ValueError, TypeError):
            pass
    genera = request.get('genera', {})
    if 'organisms' in request:
        genera = {}
        for slot, entry in request['organisms'].items():
            label = extract_genus_label(entry or {})
            genera[slot] = label.lower() if label else None
    g = [genera.get(slot) for slot in ('org1', 'org2', 'org3')]
    genus_condition = (g[0] != g[1] and g[0] != g[2] and g[1] == g[2]) if all(g) else None
    gap_active = thresholds['min_idt_gap_between_out_in'] is not None
    waived = gap_active and genus_condition is True
    statuses = {}
    reasons = []
    requirements = {
        'ingroup_idt_threshold': ['idt_23'],
        'outgroup_idt_threshold': ['idt_12', 'idt_13'],
        'min_idt_gap_between_out_in': ['idt_12', 'idt_13', 'idt_23'],
    }
    for key, required in requirements.items():
        threshold = thresholds[key]
        status = 'disabled'
        if threshold is not None:
            if key == 'min_idt_gap_between_out_in' and waived:
                status = 'waived'
            elif any(k not in values for k in required):
                status = 'missing_identity'
            else:
                if key == 'min_idt_gap_between_out_in':
                    compared = identity_gaps(values)
                else:
                    compared = [Decimal(str(values[k])) for k in required]
                status = 'passed' if all(v > Decimal(str(threshold)) for v in compared) else 'failed'
            if status in ('failed', 'missing_identity'):
                reasons.append(f'{key}: {status}')
        statuses[key] = status
    gaps = (None, None)
    if all(k in values for k in ('idt_12', 'idt_13', 'idt_23')):
        gaps = identity_gaps(values)
    return dict(
        thresholds=thresholds,
        identity_source=IDENTITY_SOURCE,
        genera=genera,
        **values,
        gap_12=float(gaps[0]) if gaps[0] is not None else None,
        gap_13=float(gaps[1]) if gaps[1] is not None else None,
        genus_condition=genus_condition,
        gap_waived=waived,
        **statuses,
        included=not reasons,
        reason='; '.join(reasons),
    )


if __name__ == '__main__':
    print(json.dumps(evaluate(json.loads(sys.stdin.read())), allow_nan=False))
