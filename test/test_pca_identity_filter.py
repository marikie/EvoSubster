import sys
from pathlib import Path
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src/lib'))
from pca_identity_filter import evaluate

class IdentityFilterTests(unittest.TestCase):
    def check_filter(self, options, identities=None, genera=None):
        return evaluate(dict(options=options, identities=identities if identities is not None else dict(idt_12=91.8,idt_13=90,idt_23=95), genera=genera or {}))

    def test_disabled_and_individual(self):
        self.assertTrue(self.check_filter({}, {})['included'])
        self.assertTrue(self.check_filter({'ingroup_idt_threshold':94},dict(idt_23=95))['included'])
        self.assertFalse(self.check_filter({'outgroup_idt_threshold':90})['included'])
        self.assertFalse(self.check_filter({'ingroup_idt_threshold':95})['included'])

    def test_common_and_override(self):
        result=self.check_filter(dict(idt_threshold=96,ingroup_idt_threshold=94,outgroup_idt_threshold=89))
        self.assertTrue(result['included'])
        self.assertFalse(self.check_filter(dict(idt_threshold=96))['included'])

    def test_exact_gap_and_waiver(self):
        option={'min_idt_gap_between_out_in':3.2}
        self.assertFalse(self.check_filter(option)['included'])
        for genera in [dict(org1='a',org2='a',org3='a'), dict(org1='a',org2='b',org3='c'), dict(org2='b',org3='b')]:
            self.assertFalse(self.check_filter(option,genera=genera)['included'])
        genera=dict(org1='a',org2='b',org3='b')
        self.assertEqual(self.check_filter(option,genera=genera)['min_idt_gap_between_out_in'],'waived')
        self.assertFalse(self.check_filter(dict(**option,ingroup_idt_threshold=95),genera=genera)['included'])
        self.assertTrue(self.check_filter(option,{},genera)['included'])

    def test_invalid_and_missing(self):
        for value in [-1,101,float('nan'),float('inf')]:
            with self.assertRaises(ValueError): self.check_filter({'idt_threshold':value})
        self.assertFalse(self.check_filter({'ingroup_idt_threshold':80},{})['included'])
        self.assertFalse(self.check_filter({'ingroup_idt_threshold':80},{'idt_23':float('nan')})['included'])

    def test_waiver_absolute_outgroup_and_passing_gap(self):
        genera=dict(org1='a',org2='b',org3='b')
        self.assertFalse(self.check_filter(dict(outgroup_idt_threshold=91,min_idt_gap_between_out_in=50),genera=genera)['included'])
        self.assertTrue(self.check_filter(dict(min_idt_gap_between_out_in=3.19))['included'])
        self.assertTrue(evaluate(dict(options={},identities=[]))['included'])

    def test_selector_genus_whitespace(self):
        result=evaluate(dict(options={'min_idt_gap_between_out_in':100},identities=[],organisms={
            'org1':{'raw_organism_name':' Pan_alpha '},
            'org2':{'raw_organism_name':' Homo\tbeta '},
            'org3':{'ncbi_full_name':'homo gamma'}}))
        self.assertTrue(result['included'])
        self.assertTrue(result['gap_waived'])
