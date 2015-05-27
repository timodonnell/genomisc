from collections import defaultdict

import pandas

from lxml import etree

def parse_to_dataframe(filename):
    raise NotImplementedError()
    parsed = etree.parse(filename)
    result = defaultdict(lambda: defaultdict(list))
    for item in parsed.getroot().findall("Result"):
        result[item.find("participant_id").text][item.find("aliquot_id").text].append(item)
    return result


