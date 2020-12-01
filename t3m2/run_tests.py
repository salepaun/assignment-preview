#! /bin/python

import os
import subprocess
import sys
import argparse

from os import path

def find_tests(args):
    """Find all of the test files in the assets directory."""
    test_files = []
    exclude = set()
    if args.omit:
        exclude = set(args.omit)
    
    assets_path = path.join(args.workpath, 'assets/{}/')
    for dirpath, dirnames, filenames in os.walk(assets_path.format(args.dir)):
        dirnames[:] = [d for d in dirnames if d not in exclude] 
        for filename in filenames:
            fullname = dirpath.split('/')[-1] + '/' + filename
            if filename.endswith('.xml') and fullname not in exclude:
                if args.specific is None or dirpath.split('/')[-1] in args.specific or fullname in args.specific:
                    test_files.append(os.path.join(dirpath, filename))

    return test_files


def main():
    """Collect the tests for the current milestone and run them against the oracle.

    If there are extra credit tests, they can be run independently of the other tests
    by using the --extra or -e flags. To run specific tests, one can use the --specific
    or -s flags with the name of the test directory:
    
    --specific -s
    --extra -e
    --omit -o

    Examples
    --------
    $  python3 run_tests.py t4m1
        $  python3 run_tests.py t4m1 -e
    $  python3 run_tests.py t4m1 -s SpringTests
    """
    
    work_path = "/home/codio/workspace";
    parser = argparse.ArgumentParser(description='Run Oracle tests')
    parser.add_argument('-s', '--specific', nargs='*')
    parser.add_argument('-e', '--extra', action='store_true')
    parser.add_argument('-o', '--omit', nargs='*')
    parser.add_argument('-p', '--workpath', default=work_path, type=str, help='Workspace path')
    parser.add_argument('dir', type=str, help='Name of the theme')
    args = parser.parse_args();
    
    theme = args.dir
    
    if args.extra:
        args.dir += "_extracredit" 
        
    tests = find_tests(args)
    successful_tests = 0
    failed_tests = 0

    for test in sorted(tests):
        print('Running test {}: ' .format(test), end='')
        test_proc = path.join(args.workpath, 'build/FOSSSim/FOSSSim')
        outfile = path.join(args.workpath, 'test_output.bin')
        oracle_proc = path.join(args.workpath, 'oracle/FOSSSimOracle{}')
        subprocess.check_output([test_proc, '-s', '{}' .format(test),
                                 '-d', '0', '-o', outfile])

        oracle = subprocess.check_output([oracle_proc.format(theme.upper()),
                                          '-s', '{}' .format(test), '-d', '0', '-i', outfile],
                                          universal_newlines=True)

        if 'Overall success: Passed.' in oracle:
            successful_tests += 1
            print('Passed.')
        elif 'Overall success: Failed.' in oracle:
            failed_tests += 1
            print('Failed.')

    print('------------------------------------------------')
    print('Successful Tests: {}' .format(successful_tests))
    print('Failed Tests: {}' .format(failed_tests))
    print('------------------------------------------------')


main()
