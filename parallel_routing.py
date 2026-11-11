"""Display-only rectangular detours for long adjacent highlighted segments."""
from collections import defaultdict
import math

EPS = 1e-7


def _legs(points):
    for a, b in zip(points, points[1:]):
        if math.dist(a, b) <= EPS:
            continue
        if abs(a[1] - b[1]) < EPS:
            axis = 0
        elif abs(a[0] - b[0]) < EPS:
            axis = 1
        else:
            raise ValueError('Parallel routing requires orthogonal segments')
        yield dict(a=tuple(a), b=tuple(b), axis=axis, level=a[1-axis],
                   lo=min(a[axis], b[axis]), hi=max(a[axis], b[axis]))


def _overlap(a, b):
    return min(a['hi'], b['hi']) - max(a['lo'], b['lo'])


def _distance(point, leg):
    axis = leg['axis']
    return math.hypot(max(leg['lo']-point[axis], 0, point[axis]-leg['hi']),
                      point[1-axis]-leg['level'])


def _detour(leg, offset, grid):
    a, b = leg['a'], leg['b']
    if abs(offset) < EPS:
        return [a, b]
    axis = leg['axis']; side = 1-axis
    stub = min(.08*grid, (leg['hi']-leg['lo'])/4)
    p, q = list(a), list(b)
    sign = 1 if b[axis] > a[axis] else -1
    p[axis] += sign*stub; q[axis] -= sign*stub
    u, v = list(p), list(q)
    u[side] += offset; v[side] += offset
    return [a, tuple(p), tuple(u), tuple(v), tuple(q), b]


def separate_parallel(paths, grid=1.0, obstacles=(), clearance=.3, min_overlap=1.0):
    """Mutate path points, keeping endpoints and fixed backbone segments intact.

    Paths provide rid, priority (2 red / 1 orange / 0 other), points and optional
    backbone. Same-reaction, same-color overlapping trunks share one offset.
    Obstacles are biological node coordinates; crossings of other lines are
    scored when choosing a side. Unplaceable groups remain unchanged and are
    recorded. Geometry, never fluxes or chemistry, is changed.
    """
    gap, threshold = clearance*grid, min_overlap*grid
    obstacles = [tuple(p) for p in obstacles]
    entries=[]; buckets=defaultdict(list)
    for index, path in enumerate(paths):
        for part, leg in enumerate(_legs(path['points'])):
            leg.update(path=index, part=part, rid=path['rid'], priority=path['priority'],
                       fixed=bool(path.get('backbone')), offset=0.0)
            entries.append(leg)
            if leg['priority'] > 0:
                buckets[(leg['rid'], leg['priority'], leg['axis'], round(leg['level']/grid, 8))].append(leg)
    groups=[]
    for key, legs in sorted(buckets.items()):
        group=[]; end=-math.inf
        for leg in sorted(legs, key=lambda x:(x['lo'],x['hi'],x['path'],x['part'])):
            if group and leg['lo'] > end+EPS:
                groups.append(group); group=[]
            group.append(leg); end=max(end,leg['hi'])
        if group: groups.append(group)
    groups.sort(key=lambda g:(not any(x['fixed'] for x in g), -g[0]['priority'],
                              g[0]['rid'],g[0]['axis'],g[0]['level'],min(x['lo'] for x in g)))
    reserved=[]; unresolved=[]; shifted=0
    def conflicts(group, offset):
        return [other for other in reserved if any(
            leg['rid'] != other['rid'] and leg['axis'] == other['axis'] and
            _overlap(leg,other) >= threshold-EPS and
            abs(leg['level']+offset-other['level']) < gap-EPS for leg in group)]
    def collision(group, offset):
        endpoints={p for leg in group for p in (leg['a'],leg['b'])}
        segments=[s for leg in group for s in _legs(_detour(leg,offset,grid))]
        return any(_distance(point, s) < .06*grid-EPS
                   for point in obstacles if point not in endpoints for s in segments)
    def crossings(group, offset):
        count=0
        for leg in group:
            for part in _legs(_detour(leg,offset,grid)):
                for other in entries:
                    if other['rid']==leg['rid'] or part['axis']==other['axis']:continue
                    h,v=(part,other) if part['axis']==0 else (other,part)
                    if h['lo']+EPS < v['level'] < h['hi']-EPS and v['lo']+EPS < h['level'] < v['hi']-EPS:
                        count+=1
        return count
    for group in groups:
        clashes=conflicts(group,0)
        offset=0.0
        if clashes:
            candidates=set()
            if not any(leg['fixed'] for leg in group):
                for leg in group:
                    for other in reserved:
                        if leg['axis']==other['axis'] and _overlap(leg,other)>=threshold-EPS:
                            for side in (-1,1):
                                candidate=other['level']+side*gap-leg['level']
                                if EPS < abs(candidate) <= 1.2*grid+EPS:candidates.add(candidate)
            viable=[value for value in candidates if not conflicts(group,value) and not collision(group,value)]
            if viable:
                # Nearest valid lane first; choose the less crowded side on ties.
                offset=min(viable,key=lambda value:(round(abs(value)/grid,8),crossings(group,value),value<0))
                shifted+=len(group)
            else:
                unresolved.append(dict(reaction_id=group[0]['rid'],axis=group[0]['axis'],
                                       level=group[0]['level'],
                                       reason='fixed backbone' if any(x['fixed'] for x in group) else 'no free lane within 1.2 grids',
                                       conflicting_reactions=sorted({x['rid'] for x in clashes})))
        for leg in group:
            leg['offset']=offset
            # Reserve the long shaft. The original endpoints remain connected by short stubs.
            reserved.append(dict(leg,level=leg['level']+offset))
    by_path=defaultdict(list)
    for leg in entries:by_path[leg['path']].append(leg)
    for index,path in enumerate(paths):
        if not any(leg['offset'] for leg in by_path[index]):continue
        points=[]
        for leg in by_path[index]:
            detour=_detour(leg,leg['offset'],grid)
            points.extend(detour if not points else detour[1:])
        path['points']=points
    return dict(clearance_grid=clearance,min_parallel_length_grid=min_overlap,
                shifted_segments=shifted,unresolved=unresolved,
                max_offset_grid=max((abs(x['offset'])/grid for x in entries),default=0),
                fixed_endpoints=True,shared_same_reaction_trunks=True)
