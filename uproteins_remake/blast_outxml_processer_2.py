"""filter out smorfs from old xml files"""


import xml.etree.ElementTree as ET
import csv

def parse_large_xml(file_path):
    try:
        with open(file_path, 'r', encoding='utf-8') as file:
            xml_content = file.read()

        root = ET.fromstring(xml_content)
        return root

    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
    except ET.ParseError as e:
        print(f"Error parsing XML: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

    return None

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

def save_keys_to_txt(dictionary, txt_file_path):
    with open(txt_file_path, mode='w', encoding='utf-8') as file:
        for key in dictionary.keys():
            file.write(f"{key}\n")

xml_file_path = '/home/adriana/uproteins/redo/phylogenetic_analysis/old_xml/tblastn/torfs_T5_2.xml'
filename = xml_file_path.replace("/home/adriana/uproteins/redo/phylogenetic_analysis/old_xml/tblastn/", "")
outfile_name = filename.replace(".xml", "")
csv_file_path = f'/home/adriana/uproteins/redo/phylogenetic_analysis/old_xml/tblastn/{outfile_name}.csv'
txt_file_path = f'/home/adriana/uproteins/redo/phylogenetic_analysis/old_xml/tblastn/{outfile_name}.txt'

root = parse_large_xml(xml_file_path)

if root is not None:
    result_dict = process_blast_output(root)
    # Print the result dictionary
    for query_def, hit_defs in result_dict.items():
        print(f"{query_def}: {hit_defs}")
    save_dict_to_csv(result_dict, csv_file_path)
    save_keys_to_txt(result_dict, txt_file_path)
    print(f"Successfully saved {filename} summary to {csv_file_path}")
    print(f"Successfully saved keys to {txt_file_path}")
else:
    print("Failed to parse the XML file.")





