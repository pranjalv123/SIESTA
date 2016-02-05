import dendropy
import argparse


taxon_namespace = None


def init(ts, opts):
    print "CALLING INIT FN"
    taxon_namespace = ts
    parser = argparse.ArgumentParser()
    print "OPTIONS", opts
def score(t):
    print "CALLING SCORE FN"
    return 0
