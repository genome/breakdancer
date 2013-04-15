#!/usr/bin/env python
"""This is a simple python script to perform an integration
test on BreakDancer. It uses some small BAMs from sequencing
of NA19238 and NA19240 and the integrationtest.py from
TGI's build-common submodule."""

from integrationtest import IntegrationTest, main
import unittest

class TestBreakDancer(IntegrationTest, unittest.TestCase):

    def test_breakdancer(self):
        expected_file = self.inputFiles("expected_output")[0]
        config_file = self.inputFiles("inv_del_bam_config")[0]
        output_file = self.tempFile("output")
        params = [
            "-o 21", " > ", output_file
        ]
        rv, err = self.execute(params)
        self.assertEqual(0, rv)
        self.assertFilesEqual(expected_file, output_file, filter_regex="#Command")

if __name__ == "__main__":
    main()
