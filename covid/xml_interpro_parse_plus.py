"""this includes pathway information"""

import xml.etree.ElementTree as ET

# Parse the XML file
tree = ET.parse('/home/adriana/Downloads/rescored_microproteins_sars.fasta.xml')
root = tree.getroot()

# Namespace dictionary
ns = {'ip': 'http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5'}

# Open a file for writing
with open('parsed_output_pathways.txt', 'w') as f:
    # Iterate through each 'protein' element
    for protein in root.findall('ip:protein', ns):
        sequence = protein.find('ip:sequence', ns).text
        xref = protein.find('ip:xref', ns)
        xref_id = xref.get('id')
        xref_name = xref.get('name')

        f.write("Sequence: {}\n".format(sequence))
        f.write("XRef ID: {}\n".format(xref_id))
        f.write("XRef Name: {}\n".format(xref_name))

        # Iterate through all matches
        for matches in protein.findall('ip:matches', ns):
            for match in matches:
                # Extract match information
                signature = match.find('ip:signature', ns)
                ac = signature.get('ac')
                desc = signature.get('desc')
                f.write("Match AC: {}\n".format(ac))
                f.write("Match Description: {}\n".format(desc))

                # Extract pathways
                pathways = signature.findall('.//ip:pathway-xref', ns)
                if pathways:
                    f.write("Pathways:\n")
                    for pathway in pathways:
                        db = pathway.get('db')
                        id = pathway.get('id')
                        name = pathway.get('name')
                        f.write("  DB: {}, ID: {}, Name: {}\n".format(db, id, name))

                # Extract additional information if needed
                locations = match.find('ip:locations', ns)
                if locations is not None:
                    for location in locations.findall('ip:mobidblite-location', ns):
                        start = location.get('start')
                        end = location.get('end')
                        f.write("Location Start: {}\n".format(start))
                        f.write("Location End: {}\n".format(end))

        f.write("\n")  # Add a line break for readability
