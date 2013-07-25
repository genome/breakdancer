configure_file(testdata.py.in testdata.py @ONLY)
def_integration_test(breakdancer-max RunBreakDancer breakdancer_test.py)
