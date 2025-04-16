import csv
import json
from warnings import warn

extractions = set()
cleaned = []
all_viruses = set()

ACRONYMS_TO_NAMES = {
    "ACLSV":"Apple chlorotic leaf spot virus",
    "ADFVd":"Apscaviroid fossulamali",
    "AGVd":"Australian grapevine viroid",
    "AHVd": "Apple hammerhead viroid-like circular RNA",
    "APLPV":"American plum line pattern virus",
    "ApMV":"Apple mosaic virus",
    "APV-1":"Asian prunus virus 1",
    "APV-2":"Asian prunus virus 2",
    "APV-3":"Asian prunus virus 3",
    "ArMV":"Arabis mosaic virus",
    "ARWV1":"Apple rubbery wood virus 1",
    "ASGV":"Apple stem grooving virus",
    "ASPV":"Foveavirus mali",
    "ASSVd":"Apple scar skin viroid",
    "BCCV-1":"Blackcurrant-associated closterovirus 1",
    "BlMaV":"Blueberry mosaic associated virus",
    "BlScV":"Blueberry scorch virus",    
    "BlShV":"Blueberry shock virus",
    "BLV":"Blueberry latent virus",
    "CGRMV":"Cherry green ring mottle virus",
    "CiVA":"Citrus virus A",
    "CMLV":"Cherry mottle leaf virus",
    "CNRMV":"Cherry necrotic rusty mottle virus",
    "CRLV":"Cherry rasp leaf virus",
    "CTLaV":"Cherry twisted leaf associated virus",
    "CVA":"Cherry virus A",
    "DdV-2":"Discula destructiva virus 2",
    "GAMaV":"Grapevine asteroid mosaic associated virus",
    "GEEV":"Grapevine endophyte endornavirus",
    "GFkV": "Grapevine fleck virus",
    "GFLV": "Grapevine fanleaf virus",
    "GLRaV1": "Grapevine leafroll-associated virus 1",
    "GLRaV2": "Grapevine leafroll-associated virus 2",
    "GLRaV3":"Grapevine leafroll-associated virus 3",
    "GLRaV4":"Grapevine leafroll-associated virus 4",
    "GLRaV10":"Grapevine leafroll-associated virus 10",
    "GPGV": "Grapevine pinot gris virus",
    "Grapevine enamovirus 1":"Grapevine enamovirus-1",
    "GRBaV":"Grapevine red-blotch associated virus",
    "GRGV": "Grapevine Red Globe Virus",
    "GRSPaV":"Grapevine rupestris stem pitting-associated virus",
    "GRVFV": "Grapevine rupestris vein feathering virus",
    "GSV":"Grapevine satellite virus",
    "GsyV-1":"Grapevine syrah virus 1",
    "GVA": "Grapevine virus A",
    "GVB":"Grapevine virus B",
    "GVBaV":"Gooseberry vein banding associated virus",
    "GVC":"Grapevine virus C",
    "GVD":"Grapevine virus D",
    "GVE":"Grapevine virus E",
    "GVF":"Grapevine virus F",
    "GVK":"Grapevine virus K",
    "GVS":"Grapevine virus S",
    "GVT":"Grapevine virus T",
    "GYSVd": "Grapevine yellow speckle viroid",
    "GYSVd 1":"Grapevine yellow speckle viroid 1",
    "HdRSV":"Hydrangea ringspot virus",
    "HSVd":"Hop stunt viroid",
    "LChV2":"Little cherry virus 2",
    "NSPaV":"Nectarine stem pitting-associated virus",
    "PBCVd":"Pear blister canker viroid",
    "PBNSPV":"Plum bark necrosis stem pitting-associated virus",
    "PcMV":"Peach mosaic virus",
    "PCMV":"Peach chlorotic mottle virus",
    "PDV": "Prune dwarf virus",
    "PLMVd":"Peach latent mosaic viroid",
    "PNRSV":"Prunus necrotic ringspot virus",
    "PpCV":"Pyrus pyrifolia cryptic virus",
    "PPV":"Plum pox virus",
    "PvEV-1":"Phaseolus vulgaris endornavirus 1",
    "PVF":"Prunus virus F",
    "RBDV":"Raspberry bushy dwarf virus",
    "RpRSV":"Raspberry ringspot virus",
    "TBRV":"Tomato black ring virus",
    "TomRSV":"Tomato ringspot virus",
}

ACRONYM_CORRECTIONS = {
    "AHvd": "AHVd",
    "Apv-1 ": "APV-1",
    "APV2": "APV-2",
    "ARMV": "ArMV",
    "`GFkV": "GFkV",
    "GFKV": "GFkV",
    "GLRaV2 ": "GLRaV2",
    "GLRAV2": "GLRaV2",
    "GPRGV":"GPGV",
    "grapevine enamovirus 1":"Grapevine enamovirus 1",
    "GRPSaV":"GRSPaV",
    "GRPSPaV":"GRSPaV",
    "GSyV-1": "GsyV-1",
    "Gsy-1":"GsyV-1",
    "GYSVd ": "GYSVd",
    "GYSVd1": "GYSVd 1",
    "HSvd": "HSVd",
    "PvEV1": "PvEV-1",
    "tomRSV": "TomRSV",
}

ACRONYM_EXCLUSIONS = {    
    "GV_TF1",
    "LV_TF1",
    
    "TyV_TF1",
}

ACRONYM_CONTAINS_EXCLUSIONS = {
    "*",
    "control",
    "dillution",
    "diluted",
    "dilution",

}

ALLOWED_MISSING_NAMES = {
    "Discula destructiva virus 2"
}


def spread_samples(sample_name:str):
    sample_name = sample_name.replace(".fastq.gz", "")

    if sample_name.endswith("1/2"):
        return [sample_name.replace("1/2", suffix) for suffix in ("1", "2")]
    
    if sample_name.endswith("A/B/C"):
        return [sample_name.replace("A/B/C", suffix) for suffix in ("A", "B", "C")]
    
    if sample_name.startswith("RZC-10"):
        return ["RZC-10" + suffix for suffix in "2/4/6/8/10/12/14/16/18/20/22/24".split("/")]        
    
    return [sample_name]

samples = []

with open("input/cleansed.csv") as f:
    reader = csv.reader(f)
    next(reader)

    for row in reader:
        name, extraction, *acronyms = row

        had_one_non_empty = False

        for acronym in acronyms:
            if not acronym or acronym in ACRONYM_EXCLUSIONS or any(part in acronym for part in ACRONYM_CONTAINS_EXCLUSIONS):
                continue

            had_one_non_empty = True
            
            try:
                virus_name = ACRONYMS_TO_NAMES[ACRONYM_CORRECTIONS.get(acronym, acronym)]
            except KeyError:
                raise ValueError(f"Unknown acronym: '{acronym}' at line {reader.line_num}")

            samples.append((name, extraction, virus_name))

        if not had_one_non_empty:
            warn(f"Sample '{name}' has no viruses at line {reader.line_num}")



with open("input/reference.json") as f:
    reference_names_to_ids = {otu["name"]: otu["_id"]
        for otu in json.load(f)["otus"]
    }

rows = []

for raw_sample_name, extraction, virus_name in samples:
    if virus_name in ALLOWED_MISSING_NAMES:
        continue

    try:
        virus_id = reference_names_to_ids[virus_name]
    except KeyError:        
        raise ValueError(f"Unknown virus name: '{name}'")
    
    for sample_name in spread_samples(raw_sample_name):    
        rows.append((sample_name, extraction, virus_name,  virus_id))


with open("input/cleansed_sanitized.csv", "w") as f:
    writer = csv.writer(f, quoting=csv.QUOTE_ALL)
    writer.writerow(("sample_name", "extraction", "virus_name", "virus_id"))

    for row in rows:
        writer.writerow(row)
