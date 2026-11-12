"""New routing rules are checked against actual geometry, not stored counts."""
import copy
import math
import unittest
from adaptive_routing import Grid,route_paths,validate_paths,segment_clearance,integer,incident_exempt
from adaptive_geometry import compact_tca
from tca_geometry import configure_plan
from sbml2map import parser,BacterialStrategy
from model_input import load_model
from pathlib import Path


def path(rid,source,target,a,b,priority):
    return dict(rid=rid,mid=target,source=source,target=target,points=[a,b],priority=priority,backbone=False)


class AdaptiveRoutingTests(unittest.TestCase):
    def test_long_uncrossed_route_allows_audited_red_crossing(self):
        paths=[path('RED','r1','r2',(1,2),(3,2),2),path('ORANGE','o1','o2',(2,1),(2,3),1)]
        xy={k:p for route in paths for k,p in [(route['source'],route['points'][0]),(route['target'],route['points'][-1])]}
        before=copy.deepcopy(xy);exceptions=route_paths(paths,xy,[-2,6,-2,6],None)
        self.assertEqual(len(exceptions),1)
        self.assertEqual(exceptions[0]['crossed_length_grid'],2.)
        self.assertGreaterEqual(exceptions[0]['uncrossed_length_grid'],6.)
        report=validate_paths(paths,xy,[-2,6,-2,6],None,exceptions)
        self.assertEqual(report['metrics']['red_orange'],1)
        self.assertEqual(xy,before)
        with self.assertRaisesRegex(ValueError,'Unaudited'):validate_paths(paths,xy,[-2,6,-2,6],None,[])

    def test_orange_branches_do_not_share_a_positive_length_trunk(self):
        paths=[path('R','r','a',(2,2),(4,2),1),path('R','r','b',(2,2),(4,3),1)]
        xy={'r':(2,2),'a':(4,2),'b':(4,3)}
        exceptions=route_paths(paths,xy,[0,6,0,6],None)
        self.assertEqual(validate_paths(paths,xy,[0,6,0,6],None,exceptions)['orange_overlap_grid'],0)

    def test_circle_clearance_checks_segment_interior(self):
        region=dict(shape='circle',center=[4,4],radius=2,node_ids=[])
        self.assertEqual(segment_clearance([0,4],[8,4],region),0.)
        self.assertAlmostEqual(segment_clearance([0,1],[8,1],region),1.)
        p=path('O','a','b',(1,2),(7,2),1)
        with self.assertRaisesRegex(ValueError,'clearance'):validate_paths([p],{'a':(1,2),'b':(7,2)},[0,8,0,8],region,[])

    def test_non_grid_lane_is_rejected(self):
        p=path('O','a','b',(1,.15),(4,.15),1)
        with self.assertRaisesRegex(ValueError,'Non-grid'):validate_paths([p],{'a':(1,.15),'b':(4,.15)},[0,6,0,6],None,[])

    def test_short_uncrossed_route_keeps_red_crossing_forbidden(self):
        paths=[path('R','a','b',(2,2),(3,2),2),path('O','c','d',(1,1),(4,3),1)]
        xy={k:p for route in paths for k,p in [(route['source'],route['points'][0]),(route['target'],route['points'][-1])]}
        exceptions=route_paths(paths,xy,[0,6,0,6],None)
        self.assertFalse(exceptions)
        self.assertEqual(validate_paths(paths,xy,[0,6,0,6],None,exceptions)['metrics']['red_orange'],0)

    def test_thin_corridor_allows_two_but_rejects_three_lines(self):
        paths=[path(str(i),'a'+str(i),'b'+str(i),(x,1),(x,4),1) for i,x in enumerate([2.,2.25,2.75])]
        for p in paths[1:]:p['auxiliary_lane']=True
        xy={k:p for route in paths for k,p in [(route['source'],route['points'][0]),(route['target'],route['points'][-1])]}
        self.assertEqual(validate_paths(paths[:2],xy,[0,6,0,6],None,[])['maximum_thin_corridor_load'],2)
        with self.assertRaisesRegex(ValueError,'Three or more'):validate_paths(paths,xy,[0,6,0,6],None,[])

    def test_shared_connection_is_confined_to_node_glyph(self):
        first=path('A','r','a',(2,2),(4,3),1);second=path('B','r','b',(2,2),(4,1),1)
        first.update(points=[(2,2),(2.25,2),(2.25,3),(4,3)],auxiliary_lane=True)
        second.update(points=[(2,2),(2.25,2),(2.25,1),(4,1)],auxiliary_lane=True)
        xy={'r':(2,2),'a':(4,3),'b':(4,1)}
        report=validate_paths([first,second],xy,[0,6,0,6],None,[])
        self.assertEqual(report['orange_overlap_grid'],0.)
        self.assertEqual(report['raw_orange_overlap_grid'],.25)
        self.assertEqual(report['node_port_nodes'],['r'])
        first['points']=[(2,2),(3,2),(3,3),(4,3)]
        second['points']=[(2,2),(3,2),(3,1),(4,1)]
        with self.assertRaisesRegex(ValueError,'Visible orange'):validate_paths([first,second],xy,[0,6,0,6],None,[])

    def test_tca_attachment_cannot_reenter_opposite_clearance_band(self):
        region=dict(shape='circle',center=[0,0],radius=2,node_ids=['ring'])
        route=path('R','ring','other',(2,0),(0,0),1);xy={'ring':(2,0),'other':(0,0)}
        self.assertTrue(incident_exempt((2,0),(0,0),route,xy,region))
        self.assertFalse(incident_exempt((2,0),(-2,0),route,xy,region))

    def test_compact_circle_anchors_have_one_integer_coordinate(self):
        base=Path(__file__).parent;model=load_model(base/'e_coli_core.xml' if (base/'e_coli_core.xml').exists() else base/'models/e_coli_core.xml');strategy=BacterialStrategy(model)
        keys={k for route in strategy.routes for mid in route['metabolite_ids'] for k in ['met:'+mid]}
        keys.update('rxn:'+r for route in strategy.routes for r in route['reaction_ids'])
        fixed,b=strategy.plan(strategy.routes,keys);fixed,b=configure_plan(fixed,b,'circle','allow');old=b['tca_geometry']['radius']
        fixed,b=compact_tca(fixed,b);region=b['tca_geometry'];self.assertAlmostEqual(region['radius']/old,2/3,delta=.005)
        for k in region['node_ids']:
            self.assertTrue(any(integer(v) for v in fixed[k]))
            self.assertAlmostEqual(math.dist(fixed[k],region['center_offset']),region['radius'])
        self.assertEqual(len(set(tuple(fixed[k]) for k in region['node_ids'])),18)
        self.assertEqual(parser().parse_args([]).crossing_policy,'adaptive-grid')
