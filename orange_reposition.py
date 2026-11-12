"""Bounded local relocation to remove positive-length orange overlaps.

Point crossings are not penalized. Shared trunks of one reaction are retained.
Every accepted candidate is fully routed and checked against the frozen red
geometry, exchange corridors, TCA exclusion and active exchange boundary.
"""
import copy
import math
from collections import defaultdict
from crossing_routing import EPS, orange_length, plan_routes
from exchange_alignment import validate_alignment
from exchange_corridors import validate_corridors
from tca_geometry import validate_exclusion


def orange_overlaps(paths):
    """Union overlap length for each distinct reaction pair; no point crossings."""
    lines = defaultdict(lambda: defaultdict(list))
    for p in paths:
        if p['priority'] != 1:
            continue
        for a, b in zip(p['points'], p['points'][1:]):
            if math.dist(a, b) < EPS:
                continue
            axis = 0 if abs(a[1]-b[1]) < EPS else 1
            lines[axis, round(a[1-axis], 8)][p['rid']].append(sorted((a[axis], b[axis])))
    overlaps = []
    for (axis, level), reactions in sorted(lines.items()):
        ids = sorted(reactions)
        for i, left in enumerate(ids):
            for right in ids[i+1:]:
                intervals = []
                for a, b in reactions[left]:
                    for c, d in reactions[right]:
                        low, high = max(a, c), min(b, d)
                        if high-low > EPS:
                            intervals.append((low, high))
                merged = []
                for low, high in sorted(intervals):
                    if merged and low <= merged[-1][1]+EPS:
                        merged[-1][1] = max(merged[-1][1], high)
                    else:
                        merged.append([low, high])
                for low, high in merged:
                    a, b = [0., 0.], [0., 0.]
                    a[axis], b[axis] = low, high
                    a[1-axis] = b[1-axis] = level
                    overlaps.append(dict(reactions=[left, right], points=[a, b], length_grid=high-low))
    return dict(length_grid=sum(p['length_grid'] for p in overlaps), segments=overlaps)


def relocation_candidates(nodes, xy, backbones, paths, radius=3):
    """Move affected ordinary endpoints; pair corridors retain their boundary axis."""
    overlaps = orange_overlaps(paths)
    affected = {r for s in overlaps['segments'] for r in s['reactions']}
    fixed = set(backbones.get('fixed_node_offsets', {}))
    fixed.update(k for p in backbones.get('red_first', {}).get('routes', []) for k in (p['source'], p['target']))
    groups = backbones.get('exchange_corridors', {}).get('groups', [])
    pairs = {g['external']: g for g in groups if g['kind'] == 'pair'}
    locked = {k for g in groups for k in g['node_ids']} - set(pairs)
    keys = {k for p in paths if p['priority'] == 1 and p['rid'] in affected for k in (p['source'], p['target'])}
    boundary = backbones['routing_boundary']
    region = backbones.get('tca_exclusion', {}).get('region')
    out = []
    for key in sorted(keys, key=lambda k: (k not in pairs, k)):
        if key in fixed | locked or nodes[key]['kind'] == 'currency':
            continue
        base = xy[key]
        group = pairs.get(key)
        axis = None
        if group:
            a, b = group['corridor']
            axis = 0 if abs(a[1]-b[1]) < EPS else 1
        offsets = [(dx, dy) for dx in range(-radius, radius+1) for dy in range(-radius, radius+1)
                   if 0 < abs(dx)+abs(dy) <= radius and (axis is None or (dy == 0 if axis == 0 else dx == 0))]
        # Upward candidates precede downward ones in ties (Cartesian y grows up).
        offsets.sort(key=lambda d: (abs(d[0])+abs(d[1]), -d[1], d[0]))
        for dx, dy in offsets:
            point = [base[0]+dx, base[1]+dy]
            if not (boundary[0]+.16 < point[0] < boundary[1]-.16 and boundary[2]+.16 < point[1] < boundary[3]-.16):
                continue
            if any(math.dist(point, p) < .16-EPS for k, p in xy.items() if k != key):
                continue
            if group:
                a = xy[group['exchange']]
                # Keep the same face and an inward-only straight pair.
                direction = base[axis]-a[axis]
                if (point[axis]-a[axis])*direction <= 0:
                    continue
            test = dict(xy);test[key] = point
            try:
                validate_exclusion(nodes, test, region)
            except (ValueError, RuntimeError):
                continue
            out.append((key, point))
    return out


def reposition_orange(nodes, edges, xy, backbones, major, groups, cols, rows,
                      attempts=4, max_candidates=24, progress=print):
    backbones.setdefault('routing_boundary', [0, cols-1, 0, rows-1])
    plan = backbones.get('crossing_routing')
    if plan is None:
        plan = plan_routes(nodes, edges, xy, backbones, major, attempts=attempts, progress=progress)
    initial = orange_overlaps(plan['routes'])
    report = dict(enabled=True, before_overlap_grid=initial['length_grid'], after_overlap_grid=initial['length_grid'],
                  before_orange_length_grid=orange_length(plan['routes']), after_orange_length_grid=orange_length(plan['routes']),
                  attempted_candidates=0, accepted_moves=[], rejected_candidates=[],
                  objective_order=['distinct_reaction_overlap_length', 'orange_length', 'displacement'],
                  point_crossings_penalized=False, shared_same_reaction_trunks_penalized=False,
                  global_optimality_proven=False, candidate_limit=max_candidates)
    tried = set()
    moved = set()
    while report['attempted_candidates'] < max_candidates:
        current_overlap = orange_overlaps(plan['routes'])['length_grid']
        if current_overlap <= EPS:
            break
        candidates = relocation_candidates(nodes, xy, backbones, plan['routes'])
        best = None
        for key, point in candidates:
            if key in moved:
                continue
            # Finish the best relocation for an affected external pair before
            # considering unrelated reactions with a larger overlap total.
            if best is not None and key != best[1]:
                break
            signature = (plan['signature'], key, tuple(point))
            if signature in tried:
                continue
            tried.add(signature)
            trial_xy = copy.deepcopy(xy);trial = copy.deepcopy(backbones)
            old = trial_xy[key];trial_xy[key] = point
            # Currency copies follow their parent, retaining their attachment distance.
            for child, node in nodes.items():
                if node['kind'] == 'currency' and node['parent_id'] == key:
                    trial_xy[child] = [xy[child][i]+point[i]-old[i] for i in (0, 1)]
            for g in trial.get('exchange_corridors', {}).get('groups', []):
                if key == g['external'] and g['kind'] == 'pair':
                    g['corridor'] = [list(trial_xy[g['exchange']]), list(point)]
            report['attempted_candidates'] += 1
            try:
                boundary = trial['routing_boundary']
                if any(not (boundary[0] <= p[0] <= boundary[1] and boundary[2] <= p[1] <= boundary[3]) for p in trial_xy.values()):
                    raise ValueError('Node outside exchange boundary')
                validate_exclusion(nodes, trial_xy, trial.get('tca_exclusion', {}).get('region'))
                validate_alignment(groups, trial_xy, cols, rows, boundary)
                candidate = plan_routes(nodes, edges, trial_xy, trial, major, attempts=attempts, progress=lambda *a, **k: None)
                validate_corridors(candidate['routes'], trial.get('exchange_corridors', {'groups': []}))
                overlap = orange_overlaps(candidate['routes'])['length_grid']
                score = (round(overlap, 7), round(orange_length(candidate['routes']), 7), math.dist(old, point))
                if overlap < current_overlap-EPS and (best is None or score < best[0]):
                    best = (score, key, list(old), list(point), trial_xy, trial, candidate)
            except (ValueError, RuntimeError) as error:
                report['rejected_candidates'].append(dict(node=key, position=point, reason=str(error)[:250]))
            if report['attempted_candidates'] >= max_candidates:
                break
        if best is None:
            break
        score, key, before, after, trial_xy, trial, plan = best
        xy.clear();xy.update(trial_xy);backbones.clear();backbones.update(trial)
        moved.add(key)
        report['accepted_moves'].append(dict(node=key, before=before, after=after, overlap_grid=score[0], orange_length_grid=score[1]))
        progress(f'Orange overlap relocation: {key} {before} -> {after}; overlap {current_overlap:.3f} -> {score[0]:.3f}', flush=True)
    backbones['crossing_routing'] = plan
    report.update(after_overlap_grid=orange_overlaps(plan['routes'])['length_grid'], after_orange_length_grid=orange_length(plan['routes']),
                  remaining_overlaps=orange_overlaps(plan['routes'])['segments'], red_routes_preserved=plan['report'].get('red_routes_preserved', False))
    backbones['orange_reposition'] = report
    return report
