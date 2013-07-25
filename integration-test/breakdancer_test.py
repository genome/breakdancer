#!/usr/bin/env python
"""This is a simple python script to perform an integration
test on BreakDancer. It uses some small BAMs from sequencing
of NA19238 and NA19240 and the integrationtest.py from
TGI's build-common submodule."""

import os
print "I AM IN", os.getcwd()
from integrationtest import IntegrationTest, main
from testdata import TEST_DATA_DIRECTORY
import unittest
import subprocess
import os

class TestBreakDancer(IntegrationTest, unittest.TestCase):

    def setUp(self):
        IntegrationTest.setUp(self)
        self.data_dir = TEST_DATA_DIRECTORY
        self.orig_path = os.path.realpath(os.getcwd())
        self.exe_path = os.path.realpath(self.exe_path)
        os.chdir(self.data_dir)

    def tearDown(self):
        IntegrationTest.tearDown(self)
        os.chdir(self.orig_path)

    def test_breakdancer_cn_per_lib(self):
        expected_file = "expected_output.cn_per_lib"
        config_file = "inv_del_bam_config"
        output_file = self.tempFile("output")
        cmdline = " ".join([self.exe_path, '-a', '-o', '21', config_file, '>', output_file])
        print "Executing", cmdline
        print "CWD", os.getcwd()
        #params = [ "-o 21", " > ", output_file ]
        #rv, err = self.execute_through_shell(params)
        rv = subprocess.call(cmdline, shell=True)
        print "Return value:", rv
        self.assertEqual(0, rv)
        self.assertFilesEqual(expected_file, output_file, filter_regex="#Command|#Software")

    def test_breakdancer_cn_per_lib_af(self):
        expected_file = "expected_output.cn_per_lib.af"
        config_file = "inv_del_bam_config"
        output_file = self.tempFile("output")
        cmdline = " ".join([self.exe_path, '-a', '-h', '-o', '21', config_file, '>', output_file])
        print "Executing", cmdline
        print "CWD", os.getcwd()
        #params = [ "-o 21", " > ", output_file ]
        #rv, err = self.execute_through_shell(params)
        rv = subprocess.call(cmdline, shell=True)
        print "Return value:", rv
        self.assertEqual(0, rv)
        self.assertFilesEqual(expected_file, output_file, filter_regex="#Command|#Software")

    def test_breakdancer_af(self):
        expected_file = "expected_output.af"
        config_file = "inv_del_bam_config"
        output_file = self.tempFile("output")
        cmdline = " ".join([self.exe_path, '-h', '-o', '21', config_file, '>', output_file])
        print "Executing", cmdline
        print "CWD", os.getcwd()
        #params = [ "-o 21", " > ", output_file ]
        #rv, err = self.execute_through_shell(params)
        rv = subprocess.call(cmdline, shell=True)
        print "Return value:", rv
        self.assertEqual(0, rv)
        self.assertFilesEqual(expected_file, output_file, filter_regex="#Command|#Software")

    def test_breakdancer(self):
        expected_file = "expected_output"
        config_file = "inv_del_bam_config"
        output_file = self.tempFile("output")
        cmdline = " ".join([self.exe_path, '-o', '21', config_file, '>', output_file])
        print "Executing", cmdline
        print "CWD", os.getcwd()
        #params = [ "-o 21", " > ", output_file ]
        #rv, err = self.execute_through_shell(params)
        rv = subprocess.call(cmdline, shell=True)
        print "Return value:", rv
        self.assertEqual(0, rv)
        self.assertFilesEqual(expected_file, output_file, filter_regex="#Command|#Software")

    def test_breakdancer_all_seqs(self):
        expected_file = "expected_output"
        config_file = "inv_del_bam_config"
        output_file = self.tempFile("output")
        cmdline = " ".join([self.exe_path, config_file, '>', output_file])
        print "Executing", cmdline
        print "CWD", os.getcwd()
        #params = [ "-o 21", " > ", output_file ]
        #rv, err = self.execute_through_shell(params)
        rv = subprocess.call(cmdline, shell=True)
        print "Return value:", rv
        self.assertEqual(0, rv)
        self.assertFilesEqual(expected_file, output_file, filter_regex="#Command|#Software")

    def test_breakdancer_all_seqs_af(self):
        expected_file = "expected_output.af"
        config_file = "inv_del_bam_config"
        output_file = self.tempFile("output")
        cmdline = " ".join([self.exe_path, '-h', config_file, '>', output_file])
        print "Executing", cmdline
        print "CWD", os.getcwd()
        #params = [ "-o 21", " > ", output_file ]
        #rv, err = self.execute_through_shell(params)
        rv = subprocess.call(cmdline, shell=True)
        print "Return value:", rv
        self.assertEqual(0, rv)
        self.assertFilesEqual(expected_file, output_file, filter_regex="#Command|#Software")

    def test_breakdancer_bed_dump(self):
        expected_file = "expected.bed"
        config_file = "inv_del_bam_config"
        output_file = self.tempFile("output")
        output_bed = self.tempFile("output")
        cmdline = " ".join([self.exe_path, '-g', output_bed, config_file, '>', output_file])
        print "Executing", cmdline
        print "CWD", os.getcwd()
        rv = subprocess.call(cmdline, shell=True)
        print "Return value:", rv
        self.assertEqual(0, rv)
        self.assertFilesEqual(expected_file, output_bed, filter_regex="#Command|#Software")

    def test_breakdancer_fastq_dump(self):
        expected_files = ["expected_output", 
                "expected.H_IJ-NA19238-NA19238-extlibs.1.fastq", 
                "expected.H_IJ-NA19238-NA19238-extlibs.2.fastq", 
                "expected.H_IJ-NA19240-NA19240-extlibs.1.fastq", 
                "expected.H_IJ-NA19240-NA19240-extlibs.2.fastq", ]
        config_file = "inv_del_bam_config"
        output_file = self.tempFile("output")
        fastq_prefix = os.path.join(self.tmp_dir,"actual")
        actual_files = [ output_file, 
                fastq_prefix + ".H_IJ-NA19238-NA19238-extlibs.1.fastq", 
                fastq_prefix + ".H_IJ-NA19238-NA19238-extlibs.2.fastq", 
                fastq_prefix + ".H_IJ-NA19240-NA19240-extlibs.1.fastq", 
                fastq_prefix + ".H_IJ-NA19240-NA19240-extlibs.2.fastq", ]
        cmdline = " ".join([self.exe_path, '-o', '21', '-d', fastq_prefix, config_file, '>', output_file])
        print "Executing", cmdline
        print "CWD", os.getcwd()
        rv = subprocess.call(cmdline, shell=True)
        print "Return value:", rv
        self.assertEqual(0, rv)
        for expected, actual in zip(expected_files, actual_files):
            self.assertFilesEqual(expected, actual, filter_regex="#Command|#Software")

if __name__ == "__main__":
    main()
