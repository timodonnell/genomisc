'''
Given a directory of TCGA datasets, make a sefara collection.

Reads from the file specfied by the -i flag (or stdin if not specified) a list
of paths to TCGA BAM files. Writes out to stdout or the file specified in
--out.

Example:

ssh demeter ls /demeter/scratch/datasets/tcga/*/TCGA-*.bam | \
    genomisc-tcga-make-sefara-collection --path-field-name demeter_nfs_path

'''
from __future__ import print_function

import argparse
import glob
import os
import pkgutil
import re
import collections
import sys

import pandas

import sefara

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO  # py 3


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out",
    help="File to write. Default: stdout")
parser.add_argument("--path-field-name", default="path",
    help="Field name to use for paths.")
parser.add_argument('-i', '--input', type=argparse.FileType('r'), default='-')
parser.add_argument("--format", choices=("python", "json"))
parser.add_argument("--path-relative", default="/",
    help="Convert paths to be relative to the given path (i.e. remove this "
    "prefix from the paths)")
parser.add_argument("--path-prepend", default="/",
    help="Prepend the given string to the paths.")

filename_pattern = re.compile(
    r'((((TCGA)-(\w{2})-(\w{4})-(\w{2})(\w)-(\w{2})(\w)-(\w{4})-(\w{2}))_?([\w-]*)?)\.bam)')
filename_pattern_fields = [
    "name",
    "barcode",
    "project",
    "tissue_source_site_code",
    "participant",
    "sample_code",
    "vial",
    "portion",
    "analyte_code",
    "plate",
    "center_code",
    "extra",
]

package_data_csv_cache = {}
def load_package_data_csv(filename):
    if filename not in package_data_csv_cache:
        raw = pkgutil.get_data("genomisc", filename)
        fd_like = StringIO(raw)
        result = pandas.read_csv(
            fd_like,
            dtype=str)
        result.index = result[result.columns[0]]
        package_data_csv_cache[filename] = result
    return package_data_csv_cache[filename]

def fields_from_filename(filename):
    match = filename_pattern.match(filename)
    if not match:
        raise ValueError("Couldn't parse filename: %s" % filename)
    result = collections.OrderedDict(
        zip(filename_pattern_fields, match.groups()[1:]))

    def load_df(name):
        return load_package_data_csv("data/tcga-code-tables/%s.csv" % name)

    # Expand 'tissue source site' field.
    df = load_df("tissueSourceSite")
    result["tissue_source_site"] = (
        df.ix[result['tissue_source_site_code']]["Source Site"])
    result["study_name"] = (
        df.ix[result['tissue_source_site_code']]["Study Name"])

    # Expand 'sample' field.
    df = load_df("sampleType")
    result["sample"] = (
        df.ix[result['sample_code']]["Definition"])

    # Expand 'analyte' field.
    df = load_df("portionAnalyte")
    result["analyte"] = (
        df.ix[result['analyte_code']]["Definition"])

    # Expand 'center' field.
    df = load_df("centerCode")
    result["center_display_name"] = (
        df.ix[result['center_code']]["Display Name"])
    result["center_short_name"] = (
        df.ix[result['center_code']]["Short Name"])

    # Parse extra field
    extra = result['extra'].lower()
    if 'illumina' in extra:
        result["platform"] = "illumina"
    elif 'solid' in extra:
        result["platform"] = "solid"
    else:
        result["platform"] = ""

    return result
    
def run():
    args = parser.parse_args()
    paths = [
        x.strip() for x in args.input.readlines() if x.strip().endswith(".bam")
    ]
    print("Found %d BAM files." % len(paths), file=sys.stderr)

    resources = []
    for path in paths:
        filename = os.path.basename(path)
        fields = fields_from_filename(filename)
        fields["name"] = "bam_%s" % fields["name"]

        resource = sefara.Resource(**fields)

        # Add some conveniences
        resource[args.path_field_name] = (
            args.path_prepend + os.path.relpath(path, args.path_relative))
        resource.uuid = os.path.basename(os.path.dirname(path))
        
        # Tags
        resource.tags.add("bam")
        if resource.platform:
            resource.tags.add(resource.platform)
        if resource.participant:
            resource.tags.add("participant_%s" % resource.participant)
        
        if "mirna" in resource.extra.lower():
            resource.tags.add("mirna")
        else:
            for analyte in ("rna", "dna"):
                if analyte in resource.analyte.lower():
                    resource.tags.add(analyte)
        for sample in ("tumor", "normal"):
            if sample in resource.sample.lower():
                resource.tags.add(sample)

        resources.append(resource)

    rc = sefara.ResourceCollection(resources)
    rc.write(file=args.out, format=args.format)
    

if __name__ == '__main__':
    run()
