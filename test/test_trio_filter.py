import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src", "select"))

import trio_filter as script

SCRIPT_PATH = os.path.join(
    os.path.dirname(__file__), "..", "src", "select", "trio_filter.py"
)

# Slot 1 is the outgroup; slots 2 and 3 are the ingroup.
OUTGROUP, INGROUP_A, INGROUP_B = "org1", "org2", "org3"


def slots(genus1, genus2, genus3):
    """Build a slot_map from three genus names (slot 1 = outgroup)."""
    return {
        OUTGROUP: {"raw_organism_name": f"{genus1} alpha"},
        INGROUP_A: {"raw_organism_name": f"{genus2} beta"},
        INGROUP_B: {"raw_organism_name": f"{genus3} gamma"},
    }


def identities(idt_12, idt_13, idt_23):
    """Outgroup-ingroup identities are idt_12/idt_13; the ingroup pair is idt_23."""
    return {"idt_12": idt_12, "idt_13": idt_13, "idt_23": idt_23}


# Identities where the ingroup pair is the most similar, so the ordering rule holds.
ORDER_OK = identities("90 %", "91 %", "97 %")
# Identities where an outgroup pair beats the ingroup pair: the ordering rule fails.
ORDER_BAD = identities("98 %", "91 %", "93 %")


class TestGenusPatterns(unittest.TestCase):
    def test_pattern1_outgroup_alone_ingroup_shares_genus(self):
        patterns, issues = script.evaluate_genus_patterns(slots("Pan", "Homo", "Homo"))
        self.assertEqual(issues, [])
        self.assertTrue(patterns["pattern1"])
        self.assertFalse(patterns["pattern3"])
        self.assertFalse(patterns["pattern4"])
        self.assertFalse(patterns["pattern5"])

    def test_pattern3_all_same_genus(self):
        patterns, _ = script.evaluate_genus_patterns(slots("Homo", "Homo", "Homo"))
        self.assertTrue(patterns["pattern3"])
        self.assertFalse(patterns["pattern1"])
        self.assertFalse(patterns["pattern5"])

    def test_pattern4_all_different_genera(self):
        patterns, _ = script.evaluate_genus_patterns(slots("Pan", "Homo", "Macaca"))
        self.assertTrue(patterns["pattern4"])
        self.assertFalse(patterns["pattern1"])
        self.assertFalse(patterns["pattern5"])

    def test_pattern5_outgroup_shares_genus_with_second_ingroup(self):
        # Outgroup is congeneric with ingroup B only: taxonomy implies ((A,C),B),
        # which contradicts the assumed outgroup assignment.
        patterns, _ = script.evaluate_genus_patterns(slots("Homo", "Pan", "Homo"))
        self.assertTrue(patterns["pattern5"])
        self.assertFalse(patterns["pattern1"])

    def test_pattern5_outgroup_shares_genus_with_first_ingroup(self):
        patterns, _ = script.evaluate_genus_patterns(slots("Homo", "Homo", "Pan"))
        self.assertTrue(patterns["pattern5"])
        self.assertFalse(patterns["pattern1"])

    def test_genus_matching_is_case_insensitive(self):
        patterns, _ = script.evaluate_genus_patterns(slots("Pan", "HOMO", "homo"))
        self.assertTrue(patterns["pattern1"])

    def test_missing_slot_reports_issue(self):
        slot_map = slots("Pan", "Homo", "Homo")
        del slot_map[INGROUP_B]
        patterns, issues = script.evaluate_genus_patterns(slot_map)
        self.assertIsNone(patterns)
        self.assertTrue(issues)


class TestIdentityCondition(unittest.TestCase):
    def test_ingroup_pair_most_similar_passes(self):
        met, issues, values = script.evaluate_identity_condition(ORDER_OK)
        self.assertTrue(met)
        self.assertEqual(issues, [])
        self.assertEqual(values["idt_23"], 97.0)

    def test_outgroup_pair_more_similar_fails(self):
        met, _, _ = script.evaluate_identity_condition(ORDER_BAD)
        self.assertFalse(met)

    def test_missing_identity_is_unknown_not_false(self):
        met, issues, values = script.evaluate_identity_condition(
            identities("90 %", None, "97 %")
        )
        self.assertIsNone(met)
        self.assertIsNone(values)
        self.assertTrue(issues)

    def test_unparseable_identity_is_unknown(self):
        met, issues, _ = script.evaluate_identity_condition(
            identities("90 %", "junk", "97 %")
        )
        self.assertIsNone(met)
        self.assertTrue(issues)

    def test_overflowing_identity_is_unknown(self):
        met, issues, values = script.evaluate_identity_condition(identities("1e999", "1e999", "1e999"))
        self.assertIsNone(met)
        self.assertIsNone(values)
        self.assertTrue(issues)


class TestFilterStatus(unittest.TestCase):
    def status(self, slot_map, identity_values, threshold=80.0, **options):
        details, _ = script.evaluate_filter_status(slot_map, identity_values, threshold, **options)
        return details

    def test_separate_absolute_thresholds_are_strict_and_independent(self):
        for key, threshold, values in (
            ("ingroup_idt_threshold", 95, (85, 86, 95)),
            ("outgroup_idt_threshold", 85, (85, 86, 95)),
            ("outgroup_idt_threshold", 86, (87, 86, 95)),
        ):
            for delta, excluded in ((-0.1, True), (0, True), (0.1, False)):
                with self.subTest(key=key, values=values, delta=delta):
                    details = self.status(
                        slots("Pan", "Homo", "Homo"),
                        identities(*(str(v + delta) for v in values)),
                        **{key: threshold},
                    )
                    self.assertEqual(details["excluded"], excluded)

    def test_both_gaps_must_strictly_exceed_the_threshold(self):
        for out1, out2, gap1, gap2, excluded in (
            (90, 91, 5, 4, False), (90, 92, 5, 3, True), (92, 90, 3, 5, True),
            (90, 92.1, 5, 2.9, True), (92.1, 90, 2.9, 5, True),
        ):
            with self.subTest(out1=out1, out2=out2):
                details = self.status(
                    slots("Homo", "Homo", "Homo"),
                    identities(str(out1), str(out2), "95"),
                    min_idt_gap_between_out_in=3,
                )
                self.assertEqual(details["excluded"], excluded)
                self.assertEqual(details["idt_gap_condition"], not excluded)
                self.assertEqual(details["idt_gap_12"], gap1)
                self.assertEqual(details["idt_gap_13"], gap2)

    def test_fractional_gap_boundaries_do_not_pass_due_to_roundoff(self):
        for out1, out2, excluded in (
            ("91.8", "91.8", True),
            ("91.8", "91.7999", True),
            ("91.7999", "91.8", True),
            ("91.8001", "91.7999", True),
            ("91.7999", "91.8001", True),
            ("91.7999", "91.7999", False),
        ):
            with self.subTest(out1=out1, out2=out2):
                details = self.status(
                    slots("Homo", "Homo", "Homo"), identities(out1, out2, "95"),
                    min_idt_gap_between_out_in=3.2,
                )
                self.assertEqual(details["excluded"], excluded)
                self.assertEqual(details["idt_gap_condition"], not excluded)
                if out1 == "91.8":
                    self.assertEqual(details["idt_gap_12"], 3.2)
                if out2 == "91.8":
                    self.assertEqual(details["idt_gap_13"], 3.2)

    def test_genus_exemption_applies_to_configured_gap_only(self):
        options = dict(ingroup_idt_threshold=92, outgroup_idt_threshold=90,
                       min_idt_gap_between_out_in=10)
        details = self.status(slots("Pan", "Homo", "Homo"), ORDER_BAD, **options)
        self.assertFalse(details["excluded"])
        self.assertTrue(details["idt_gap_exempt"])
        self.assertFalse(details["idt_gap_condition"])
        for values in (identities("98", "90", "93"), identities("98", "91", "92")):
            self.assertTrue(self.status(slots("Pan", "Homo", "Homo"), values, **options)["excluded"])
        self.assertTrue(self.status(slots("Homo", "Homo", "Pan"), ORDER_OK, **options)["excluded"])

    def test_default_zero_gap_rejects_equal_identity_without_exemption(self):
        details = self.status(slots("Homo", "Homo", "Homo"), identities("95", "90", "95"))
        self.assertTrue(details["excluded"])

    def test_pattern1_retained_even_when_identity_ordering_fails(self):
        # The rule that distinguishes the thesis method: the target genus
        # configuration is always retained, without consulting the identities.
        details = self.status(slots("Pan", "Homo", "Homo"), ORDER_BAD)
        self.assertFalse(details["excluded"])
        self.assertTrue(details["genus_condition_met"])
        self.assertFalse(details["identity_condition_met"])

    def test_all_same_genus_retained_when_ordering_holds(self):
        details = self.status(slots("Homo", "Homo", "Homo"), ORDER_OK)
        self.assertFalse(details["excluded"])

    def test_all_different_genera_retained_when_ordering_holds(self):
        details = self.status(slots("Pan", "Homo", "Macaca"), ORDER_OK)
        self.assertFalse(details["excluded"])

    def test_two_vs_one_excluded_even_when_ordering_holds(self):
        # Identities look fine, but the genera contradict the outgroup assignment.
        details = self.status(slots("Homo", "Homo", "Pan"), ORDER_OK)
        self.assertTrue(details["excluded"])
        self.assertTrue(details["genus_pattern_two_vs_one"])
        self.assertTrue(details["identity_condition_met"])

    def test_failed_ordering_excludes_when_genus_is_uninformative(self):
        details = self.status(slots("Homo", "Homo", "Homo"), ORDER_BAD)
        self.assertTrue(details["excluded"])

    def test_threshold_gate_excludes_before_other_rules(self):
        # Below-threshold identities are rejected even in the target genus
        # configuration, which would otherwise be retained unconditionally.
        low = identities("70 %", "71 %", "97 %")
        details = self.status(slots("Pan", "Homo", "Homo"), low)
        self.assertTrue(details["excluded"])
        self.assertFalse(details["idt_threshold_condition"])

    def test_threshold_is_strict(self):
        at_threshold = identities("80 %", "90 %", "97 %")
        self.assertTrue(self.status(slots("Pan", "Homo", "Homo"), at_threshold)["excluded"])
        just_above = identities("80.1 %", "90 %", "97 %")
        self.assertFalse(self.status(slots("Pan", "Homo", "Homo"), just_above)["excluded"])

    def test_threshold_is_configurable(self):
        values = identities("85 %", "86 %", "97 %")
        self.assertFalse(self.status(slots("Pan", "Homo", "Homo"), values, 80.0)["excluded"])
        self.assertTrue(self.status(slots("Pan", "Homo", "Homo"), values, 90.0)["excluded"])

    def test_missing_identity_leaves_verdict_unknown_and_not_excluded(self):
        details = self.status(slots("Homo", "Homo", "Homo"), identities(None, None, None))
        self.assertIsNone(details["filter_flag"])
        self.assertFalse(details["excluded"])


class TestBatchCLI(unittest.TestCase):
    HEADER = "out_acc\tin1_acc\tin2_acc\tidt_12\tidt_13\tidt_23\tgenus_1\tgenus_2\tgenus_3"

    def run_cli(self, rows, extra_args=()):
        with tempfile.TemporaryDirectory() as tmp:
            trios = os.path.join(tmp, "trios.tsv")
            out = os.path.join(tmp, "verdict.tsv")
            with open(trios, "w", encoding="utf-8") as handle:
                handle.write(self.HEADER + "\n")
                for row in rows:
                    handle.write("\t".join(row) + "\n")
            result = subprocess.run(
                [sys.executable, SCRIPT_PATH, "--trios", trios, "--out", out, *extra_args],
                capture_output=True,
                text=True,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            with open(out, encoding="utf-8") as handle:
                lines = handle.read().strip().split("\n")
        header = lines[0].split("\t")
        return [dict(zip(header, line.split("\t"))) for line in lines[1:]]

    def test_retained_and_excluded_rows(self):
        rows = [
            # pattern1, ordering fails -> retained anyway
            ("GCA_1", "GCA_2", "GCA_3", "98", "91", "93", "Pan", "Homo", "Homo"),
            # two-vs-one, ordering holds -> excluded
            ("GCA_4", "GCA_5", "GCA_6", "90", "91", "97", "Homo", "Homo", "Pan"),
            # below threshold -> excluded
            ("GCA_7", "GCA_8", "GCA_9", "70", "71", "97", "Pan", "Homo", "Homo"),
        ]
        verdicts = self.run_cli(rows)
        self.assertEqual([v["excluded"] for v in verdicts], ["FALSE", "TRUE", "TRUE"])
        self.assertEqual(verdicts[1]["genus_pattern_two_vs_one"], "TRUE")
        self.assertEqual(verdicts[2]["idt_threshold_condition"], "FALSE")

    def test_input_columns_are_passed_through(self):
        rows = [("GCA_1", "GCA_2", "GCA_3", "90", "91", "97", "Pan", "Homo", "Homo")]
        verdict = self.run_cli(rows)[0]
        self.assertEqual(verdict["out_acc"], "GCA_1")
        self.assertEqual(verdict["genus_3"], "Homo")

    def test_threshold_option_changes_verdict(self):
        rows = [("GCA_1", "GCA_2", "GCA_3", "85", "86", "97", "Pan", "Homo", "Homo")]
        self.assertEqual(self.run_cli(rows)[0]["excluded"], "FALSE")
        strict = self.run_cli(rows, ("--idt-threshold", "90"))
        self.assertEqual(strict[0]["excluded"], "TRUE")

    def test_individual_options_override_common_threshold_in_any_order(self):
        rows = [("GCA_1", "GCA_2", "GCA_3", "85", "86", "97", "Homo", "Homo", "Homo")]
        common = ("--idt-threshold", "99")
        individual = ("--ingroup-idt-threshold", "95", "--outgroup-idt-threshold", "84",
                      "--min-idt-gap-between-out-in", "10")
        for args in (common + individual, individual + common):
            verdict = self.run_cli(rows, args)[0]
            self.assertEqual(verdict["excluded"], "FALSE")
            self.assertEqual(float(verdict["ingroup_idt_threshold_value"]), 95)
            self.assertEqual(float(verdict["outgroup_idt_threshold_value"]), 84)
            self.assertEqual(float(verdict["min_idt_gap_between_out_in_value"]), 10)
            self.assertEqual(float(verdict["idt_gap_12"]), 12)
            self.assertEqual(float(verdict["idt_gap_13"]), 11)
            self.assertEqual(verdict["idt_gap_exempt"], "FALSE")
        self.assertEqual(self.run_cli(rows, ("--min-idt-gap-between-out-in", "11"))[0]["excluded"], "TRUE")

    def test_default_thresholds_are_recorded(self):
        rows = [("GCA_1", "GCA_2", "GCA_3", "90", "91", "97", "Pan", "Homo", "Homo")]
        verdict = self.run_cli(rows)[0]
        self.assertEqual(float(verdict["ingroup_idt_threshold_value"]), 80)
        self.assertEqual(float(verdict["outgroup_idt_threshold_value"]), 80)
        self.assertEqual(float(verdict["min_idt_gap_between_out_in_value"]), 0)
        self.assertEqual(verdict["idt_gap_exempt"], "TRUE")

    def test_invalid_thresholds_fail_before_reading_input(self):
        for option in ("--idt-threshold", "--ingroup-idt-threshold", "--outgroup-idt-threshold",
                       "--min-idt-gap-between-out-in"):
            for value in ("-1", "100.1", "nan", "inf", "junk", None):
                with self.subTest(option=option, value=value):
                    result = subprocess.run(
                        [sys.executable, SCRIPT_PATH, "--trios", "/nonexistent/trios.tsv",
                         option, *([] if value is None else [value])], capture_output=True, text=True,
                    )
                    self.assertNotEqual(result.returncode, 0)
                    self.assertIn(option, result.stderr)
                    self.assertNotIn("FileNotFoundError", result.stderr)

    def test_unknown_verdict_is_written_as_na(self):
        rows = [("GCA_1", "GCA_2", "GCA_3", "", "", "", "Homo", "Homo", "Homo")]
        verdict = self.run_cli(rows)[0]
        self.assertEqual(verdict["filter_flag"], "NA")
        self.assertEqual(verdict["excluded"], "FALSE")

    def test_missing_required_column_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            trios = os.path.join(tmp, "trios.tsv")
            with open(trios, "w", encoding="utf-8") as handle:
                handle.write("out_acc\tidt_12\tidt_13\tidt_23\n")
                handle.write("GCA_1\t90\t91\t97\n")
            result = subprocess.run(
                [sys.executable, SCRIPT_PATH, "--trios", trios],
                capture_output=True,
                text=True,
            )
        self.assertEqual(result.returncode, 1)
        self.assertIn("genus_1", result.stderr)


if __name__ == "__main__":
    unittest.main()
