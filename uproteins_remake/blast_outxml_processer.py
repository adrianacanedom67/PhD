"""parse xml blast output file - gather only those that have matches and have a bit score > 50"""


import xml.etree.ElementTree as ET
import csv



def parse_large_xml(file_path):
    try:
        # Attempt to open the file
        # print("Attempting to open the file...")
        with open(file_path, 'r', encoding='utf-8') as file:
            xml_content = file.read()
            # print("Successfully read file")

        # Optionally, add a debug print to show a snippet of the file content
        # print("File content snippet:", xml_content[:1000])  # Adjust the snippet size as needed

        # Parse the XML content directly
        # print("Attempting to parse the XML content...")
        root = ET.fromstring(xml_content)
        # print("Successfully parsed XML content")
        return root

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
    except ET.ParseError as e:
        print(f"Error parsing XML: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return None

# Function to process the XML and create the dictionary
def process_blast_output(xml_root):
    result_dict = {}

    for iteration in xml_root.findall('.//Iteration'):
        query_def = iteration.find('Iteration_query-def').text
        hits = []

        for hit in iteration.findall('.//Hit'):
            hit_def = hit.find('Hit_def').text
            hsp = hit.find('.//Hsp')
            if hsp is not None:
                bit_score = float(hsp.find('Hsp_bit-score').text)
                if bit_score > 50:
                    hits.append(hit_def)

        if hits:
            result_dict[query_def] = hits

    return result_dict


def save_dict_to_csv(dictionary, csv_file_path):
    with open(csv_file_path, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        for key, values in dictionary.items():
            row = [key] + values
            writer.writerow(row)


xml_file_path = '/home/adriana/uproteins/redo/phylogenetic_analysis/blastp/incMyco_excMtb.xml'
filename = xml_file_path.replace("/home/adriana/uproteins/redo/phylogenetic_analysis/blastp/", "")
outfile_name = filename.replace(".xml", "")
csv_file_path = f'/home/adriana/uproteins/redo/phylogenetic_analysis/blastp/{outfile_name}.csv'
root = parse_large_xml(xml_file_path)

if root is not None:
    result_dict = process_blast_output(root)
    # Print the result dictionary
    for query_def, hit_defs in result_dict.items():
        print(f"{query_def}: {hit_defs}")
    save_dict_to_csv(result_dict, csv_file_path)
    print(f"succesfully saved {filename} summary to {csv_file_path}")
else:
    print("Failed to parse the XML file.")