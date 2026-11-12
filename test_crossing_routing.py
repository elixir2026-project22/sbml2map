"""Crossing constraints include touching, overlaps and circular Beziers."""
import copy
import unittest
from crossing_routing import measure,route
from tca_curves import arc_controls


def path(rid,a,b,start,end,priority=2,points=None):
    return dict(rid=rid,mid=b,source=a,target=b,priority=priority,backbone=False,points=points or [start,end])


class CrossingRoutingTests(unittest.TestCase):
    def test_unrelated_crossings_contacts_and_overlaps_are_hard(self):
        a=path('A','a','b',(0,1),(2,1))
        for p in [path('B','c','d',(1,0),(1,2)),path('B','c','d',(1,1),(1,2)),path('B','c','d',(1,1),(3,1))]:
            self.assertFalse(measure([a,p])['hard_verified'])
        p['priority']=1
        self.assertGreater(measure([a,p])['red_orange'],0)

    def test_common_biological_endpoint_and_connected_reaction_trunk(self):
        a=path('R','r','a',(0,0),(2,0))
        b=path('S','a','b',(2,0),(2,1))
        self.assertTrue(measure([a,b])['hard_verified'])
        c=path('R','r','c',(0,0),(1,1),points=[(0,0),(1,0),(1,1)])
        self.assertTrue(measure([a,c])['hard_verified'])
        # Same reaction alone cannot excuse a disconnected crossing.
        c=path('R','r','c',(0,0),(2,2),points=[(0,0),(0,-1),(3,-1),(3,1),(1,1),(1,2),(2,2)])
        a=path('R','r','a',(0,0),(2,0),points=[(0,0),(2,0)])
        self.assertTrue(measure([a,c])['hard_verified'])
        c['points']=[(0,0),(0,-1),(1,-1),(1,2),(2,2)]
        self.assertFalse(measure([a,c])['hard_verified'])

    def test_curve_interior_and_tangency_are_detected(self):
        curve=path('TCA','a','b',(1,0),(0,1));curve['curve_controls']=arc_controls((1,0),(0,1),(0,0))
        line=path('B','c','d',(-1,.5),(2,.5),priority=1)
        self.assertGreater(measure([curve,line])['red_orange'],0)
        line['points']=[(-1,1),(2,1)]
        self.assertGreater(measure([curve,line])['red_orange'],0)

    def test_red_and_orange_reroute_without_crossings(self):
        for priority in (1,2):
            paths=[path('A','a','b',(0,1),(2,1)),path('B','c','d',(1,0),(1,2),priority=priority)]
            xy={k:p for q in paths for k,p in [(q['source'],q['points'][0]),(q['target'],q['points'][-1])]}
            before=copy.deepcopy(xy)
            result=route(paths,xy,attempts=2,progress=lambda *a,**k:None)
            self.assertTrue(result['verified']);self.assertTrue(measure(paths)['hard_verified'])
            self.assertEqual(xy,before)
            if priority==1:self.assertEqual(result['metrics']['orange_orange'],0)

    def test_closed_curve_is_a_real_barrier(self):
        points=[(1,0),(0,1),(-1,0),(0,-1)];paths=[];xy={}
        for i,(a,b) in enumerate(zip(points,points[1:]+points[:1])):
            p=path('R'+str(i),'n'+str(i),'n'+str((i+1)%4),a,b);p['curve_controls']=arc_controls(a,b,(0,0));paths.append(p)
            xy[p['source']]=a;xy[p['target']]=b
        p=path('X','inside','outside',(0,0),(2,0),priority=1);paths.append(p);xy.update(inside=(0,0),outside=(2,0))
        result=route(paths,xy,attempts=0,progress=lambda *a,**k:None)
        self.assertFalse(result['verified'])

    def test_orange_crossing_is_reported_without_hard_failure(self):
        paths=[path('A','a','b',(0,1),(2,1),priority=1),path('B','c','d',(1,0),(1,2),priority=1)]
        r=measure(paths);self.assertTrue(r['hard_verified']);self.assertEqual(r['orange_orange'],1)
