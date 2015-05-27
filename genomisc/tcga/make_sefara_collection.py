'''
Given a directory of TCGA datasets, make a sefara collection.

'''

import argparse
import glob
import os
import pkgutil
import re
import collections

import pandas

import sefara

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO  # py 3

parser = argparse.ArgumentParser()
parser.add_argument("dir")
parser.add_argument("--out")
parser.add_argument("--path-relative", default="/")
parser.add_argument("--path-field-name", default="path")
parser.add_argument("--format", choices=("python", "json"))

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
    paths = glob.glob(os.path.join(args.dir, "*/TCGA-*.bam"))
    print("Found %d BAM files." % len(paths))

    resources = []
    for path in paths:
        filename = os.path.basename(path)
        fields = fields_from_filename(filename)

        # Add a few more fields.
        fields["tags"] = ["bam"]
        fields['uuid'] = os.path.basename(os.path.dirname(path))

        fields[args.path_field_name] = (
            os.path.relpath(path, args.path_relative))

        resource = sefara.Resource(**fields)
        resources.append(resource)

    rc = sefara.ResourceCollection(resources)
    rc.write(file=args.out, format=args.format)
    

if __name__ == '__main__':
    run()
