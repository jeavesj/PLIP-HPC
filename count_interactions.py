import xml.etree.ElementTree as ET
from typing import Dict

def count_interactions(xmlfile: str) -> Dict[str, int]:
    
    tree_root = ET.parse(xmlfile).getroot()

    int_types = ['hydrophobic_interaction', 'hydrogen_bond', 'water_bridge', 'salt_bridge', 'pi_stack', 'pi_cation_interaction', 'halogen_bond']
    
    count_dict = {}
    for int_type in int_types:
        count_dict[int_type] = len(tree_root.findall(f'./bindingsite/interactions/{int_type}s/{int_type}'))
    
    return count_dict
