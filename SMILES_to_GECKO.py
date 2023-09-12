"""Functions to convert to and from the GECKO structural format."""
#imports
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
import itertools
import warnings

def H_format(H_num):
    if H_num > 1:
        return f"H{H_num}"
    elif H_num == 1:
        return "H"
    else:
        return ""

def depth(obj):
    """Checks the depth of an iterable"""
    xs = (
        obj if isinstance(obj, (tuple, list, set)) else
        obj.values() if isinstance(obj, dict) else
        None
    )
    return (
        0 if xs is None else           # Not a collection.
        1 if not xs else               # Empty collection.
        1 + max(depth(x) for x in xs)  # Non-empty.
    )

def GetAtomWithIntProp(mol, IntPropValue, IntPropName="origIdx"):
    """Returns the first atom in mol with a matching Int property."""
    for a in mol.GetAtoms():
        if a.GetIntProp(IntPropName) == IntPropValue:
            return a
    return None

def IdxToIntProp(mol, Idx, IntPropName="origIdx"):
    """Converts provided indices for atoms in mol to their int property values.
    Can pass individual indices or a list/tuple with a max depth of 2."""
    if type(Idx) == int:
        return mol.GetAtomWithIdx(Idx).GetIntProp(IntPropName)
    
    elif (type(Idx) == list) and (depth(Idx) == 1):
        return [mol.GetAtomWithIdx(i).GetIntProp(IntPropName) for i in Idx]
    
    elif (type(Idx) == tuple) and (depth(Idx) == 1):
        return tuple(mol.GetAtomWithIdx(i).GetIntProp(IntPropName) for i in Idx)
    
    elif (type(Idx) == list) and (depth(Idx) == 2):
        out_lst = []
        for sublist in Idx:
            out_lst.append([mol.GetAtomWithIdx(i).GetIntProp(IntPropName) for i in sublist])
        return out_lst
    
    elif (type(Idx) == tuple) and (depth(Idx) == 2):
        out_tup_lst = []
        for subtup in Idx:
            out_tup_lst.append(tuple(mol.GetAtomWithIdx(i).GetIntProp(IntPropName) for i in subtup))
        return tuple(x for x in out_tup_lst)
    
    else:
        raise Exception(f"Unrecognised index type: {type(Idx)} with depth {depth(Idx)}.")        

def Find_Root_Frag(frags, root_orig_idx, IntPropName="origIdx"):
    """Returns a molecule that contains a given IntProp from a set of molecules"""
    for frag in frags:
        init_idxs = []
        for x in frag.GetAtoms():
            if x.HasProp("origIdx"):
                init_idxs.append(IdxToIntProp(frag, x.GetIdx(), IntPropName))
        if root_orig_idx in init_idxs:
            return frag
    raise Exception("None of the provided fragments contained the root atom.")

def SMILES_to_Mol(smiles):
    """Converts a SMILES string to an rdkit molecule with each atom labeled with
    a custom int property."""
    #read in the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    mol = rdmolops.AddHs(mol)
    
    #label all of the atoms with their initial indices to allow them to be 
    #tracked when fragmenting and editing the molecule
    for a in mol.GetAtoms():
        a.SetIntProp("origIdx", a.GetIdx())
        
    return mol

def Mol_to_GroupDicts(mol, fragment = False):
    """Takes a SMILES string and reutrns the following:
        - The ring connection positions
        - The double bond positions
        - The side-chains attached to a 'backbone' of C and O atoms
        - The functional groups each 'backbone' atom
        - The H atoms attached to each 'backbone' atom"""
       
    backbone_subunit = "[$([#6]),$([OX2]([#6])[#6])]"
    fragment_label = "[#0*]"
    
    #if this is a fragment then give the dummy atom from fragmentation an origIdx value and record it
    if fragment:
        frag_atom_curr_i = mol.GetSubstructMatch(Chem.MolFromSmarts(fragment_label))[0]
        frag_atom = mol.GetAtomWithIdx(frag_atom_curr_i)
        frag_atom.SetIntProp("origIdx", -1)
        frag_atom_i = -1
    else:
        fragment_label = ""
        
    #check for C=C double bonds in the molecule
    double_cs = [] #list to store double bonds between carbon atoms
    db_atoms = IdxToIntProp(mol,
                            mol.GetSubstructMatches(Chem.MolFromSmarts("[#6]=[#6]")))
    for a1, a2 in db_atoms:
        double_cs.append((a1,a2))
            
    #identify any atoms that are part of a ring and assign a joining point for each ring
    ring_atom_dict = {i+1 : list(ring_atoms) for i,ring_atoms in enumerate(rdmolops.GetSSSR(mol))}
    split_bonds = {}
    ring_joins = []
    for ring_no, ring_atoms in ring_atom_dict.items():
        other_ring_nos = [x for x in ring_atom_dict if x != ring_no]
        other_ring_atoms = set()
        for o_ring_no in other_ring_nos:
            for o_ring_atom in ring_atom_dict[o_ring_no]:
                other_ring_atoms.add(o_ring_atom)
        #find one of the ring atoms that is unique to this ring
        for i in ring_atoms:
            if i not in other_ring_atoms:
                unique_atom_i = i
                unique_atom = mol.GetAtomWithIdx(unique_atom_i)
        #find the ring neighbour for this unique atom
        unique_neighs = unique_atom.GetNeighbors()
        for a in unique_neighs:
            if a.GetIdx() in ring_atoms:
                #avoid splitting on double bonds
                isDouble = mol.GetBondBetweenAtoms(unique_atom_i, a.GetIdx()).GetBondType() == Chem.rdchem.BondType.DOUBLE
                if not isDouble:
                    neigh_atom = a
                    neigh_atom_i = a.GetIdx()
                
        split_bonds[ring_no] = mol.GetBondBetweenAtoms(unique_atom_i, neigh_atom_i).GetIdx()
        ring_joins.append((unique_atom.GetIntProp("origIdx"), ring_no))
        ring_joins.append((neigh_atom.GetIntProp("origIdx"), ring_no))
        
    #Split the rings open for the rest of the process
    if split_bonds:
        mol = rdmolops.FragmentOnBonds(mol, split_bonds.values(), 
                                       addDummies = False)

    #find the longest "backbone" in the moleule (C chain or ether oxygens)
    if fragment:
        pat_str = fragment_label + "~" + backbone_subunit
    else:
        pat_str = backbone_subunit
    pat = Chem.MolFromSmarts(pat_str)
    long_idxs = IdxToIntProp(mol, mol.GetSubstructMatch(pat))
    
    if long_idxs: #ensure the provided species has a carbon (or oxygen) chain
        match_bool = True
    else:
        raise Exception("Provided SMILES does not include a 'backbone' (C or O)")
    
    #iteratively increase the chain length until there is no longer a match
    while match_bool:
        pat_str += "~"+backbone_subunit #add another C or O to the backbone
        
        pat = Chem.MolFromSmarts(pat_str)
        matches = IdxToIntProp(mol, mol.GetSubstructMatches(pat))
        
        if matches: #check if there is a chain of this length
            long_idxs = matches[0]
            #TODO set any preferences for which backbone chain is selected
        else:
            match_bool = False
    
    #remove wildcards from frgmentation from the backbone list
    if fragment:
        long_idxs = tuple(x for x in long_idxs if x != frag_atom_i)
    
    #dicts and lists to store info about the backbone atoms
    backbone_hs = {a : [] for a in long_idxs}
    backbone_chains = {a : [] for a in long_idxs}
    backbone_grps = {a : [] for a in long_idxs}
    
    #go through each atom in the backbone and store the number of H atoms and attached groups
    for root_atom_i in long_idxs:
        root_atom = GetAtomWithIntProp(mol, root_atom_i)
        
        #Split each side-group off (i.e. any adjacent atom that isn't also in the backbone)
        side_grp_idxs = []
        chain_neighbours = []
        h_count = 0
        for neigh_atom in root_atom.GetNeighbors():
            neigh_atom_i = neigh_atom.GetIntProp("origIdx")
            
            if (neigh_atom.GetSymbol() == "H"): #count H atoms
                h_count += 1
            elif fragment and (neigh_atom_i == frag_atom_i): #ignore the fragmentation wildcard if present
                pass
            elif neigh_atom_i not in long_idxs:
                side_grp_idxs.append(neigh_atom_i)
            else:
                chain_neighbours.append(neigh_atom_i)
                
        #Record the number of H atoms
        backbone_hs[root_atom_i] = H_format(h_count)
        
        #split this backbone atom from the rest of the chain (if it is longer than 1 atom)
        if len(long_idxs) > 1:
            backbone_bonds = [mol.GetBondBetweenAtoms(root_atom.GetIdx(), 
                                                      GetAtomWithIntProp(mol, nb).GetIdx()).GetIdx() for nb in chain_neighbours]
            split_mol = rdmolops.FragmentOnBonds(mol, backbone_bonds, 
                                                 addDummies = True, 
                                                 dummyLabels = [(0,0)]*len(backbone_bonds))
            sep_splits = rdmolops.GetMolFrags(split_mol, asMols=True, 
                                              sanitizeFrags = False)
            
            #find the root atom fragment
            root_frag = Find_Root_Frag(sep_splits, root_atom_i)

        else:
            root_frag = mol
 
        #detect any functional groups on this atom
        #Note that many of the SMARTS strings here come from https://www.daylight.com/dayhtml_tutorials/languages/smarts/smarts_examples.html#N
        
        # detect side chains that need to be recursively run through this function
        while root_frag.GetSubstructMatch(Chem.MolFromSmarts("[#0*]"+ backbone_subunit+"~"+backbone_subunit)):
            #fragment the side-chain from the root fragment to allow the recursion
            #to work and also to allow subsequent functionality checks to ignore
            #the side-chain
            chain_match = root_frag.GetSubstructMatch(Chem.MolFromSmarts("[#0*]"+ backbone_subunit+"~"+backbone_subunit))
            chain_bond = root_frag.GetBondBetweenAtoms(chain_match[1], chain_match[2]).GetIdx()
            
            split_chain = rdmolops.FragmentOnBonds(root_frag, [chain_bond], 
                                                   addDummies = True, 
                                                   dummyLabels = [(0,0)])
            sep_split_chain = rdmolops.GetMolFrags(split_chain, asMols=True, 
                                                   sanitizeFrags = False)
            
            root_frag = Find_Root_Frag(sep_split_chain, root_atom_i)
            root_frag_i = sep_split_chain.index(root_frag)
            if root_frag_i == 0:
                side_chain = sep_split_chain[1]
            else:
                side_chain = sep_split_chain[0]
            
            side_chain_dict = Mol_to_GroupDicts(side_chain, True)
            backbone_chains[root_atom_i].append(side_chain_dict)
       
        #dictionary of SMARTS for each functional group and the GECKO group to 
        #add to the final output string
        #TODO Add R-NO2 group (be sure it doesn't interfere with R-ONO2)
        group_smarts = {"[$([#6][OX2H1])!$([#6](=O))]":"(OH)", #alcohol (avoiding carboxylic acid)
                        "[$([#6][OX2][OX2H1])!$([#6](=O))]":"(OOH)", #hydroperoxide (avoiding peracid)
                        "[$([#6][OX2][NX3+]([OX1-])(=[OX1]))!$([#6](=O))]" : "(ONO2)", #organonitrate (excluding PAN)
                        "[$([#6][O][OX1])!$([#6](=O))]":"(OO.)", #Peroxy radical (avoiding acyl peroxy)
                        "[$([#6][OX1])!$([#6](=O))]" : "(O.)", # alkoxy radical (avoiding acyl alkoxy)
                        "[CX3H1](=O)" : "O", # aldehyde
                        "[#6,#0*][CX3](=O)[#6,#0*]" : "O", #ketone
                        "[CX3](=O)[OX2H1]" : "O(OH)", #carboxylic acid
                        "[CX3](=O)[OX2][OX2H1]" : "O(OOH)", #peracid
                        "[CX3](=O)[OX2][OX1]" : "O(OO.)", #acyl peroxide
                        "[$([CX3](=O)[OX2][OX2][NX3+]([OX1-])=[OX1])]" : "O(OONO2)" #PAN
                        }
        #dictionry of smarts that should raise an error (e.g. because they are not suppported by GECKO)
        illegal_smarts = {"[CX3](=O)[OX1]" : "Acyl alkoxy radical",
                          "[F,Cl,Br,I]" : "Halogen"}
        
        for smarts, grp_name in illegal_smarts.items():
            if root_frag.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                raise Exception(f"Unsupported functional group detected: {grp_name}.")

        for smarts, gecko_grp in group_smarts.items():
            if root_frag.HasSubstructMatch(Chem.MolFromSmarts(smarts)):
                backbone_grps[root_atom_i].insert(0, gecko_grp) 
    return ring_joins, double_cs, backbone_chains, backbone_grps, backbone_hs

def GroupDict_to_GECKO(mol, ringList, dbList, chainDict, grpDict, hDict,
                       compl_atoms = set(), chain_parent = -1):
    """Converts a provided dictionary of groups attached to a 'backbone' of C 
    and O atoms into a string able to be entered into GECKO."""
    #go through each backbone atom and build up the output string
    gecko_str = ""
    backbone_atoms = grpDict.keys()
    
    compl_dbs = []
    all_db_atoms = [item for sublist in dbList for item in sublist]
    for root_atom_i in backbone_atoms:
        root_atom = GetAtomWithIntProp(mol, root_atom_i)
        
        if root_atom.GetSymbol() == "C":            
            #check if this is a double bond C for future adjustments
            is_double = root_atom_i in all_db_atoms
            if is_double:
                db_pair = [x for x in dbList if root_atom_i in x][0]
            
            #add the double bond symbol if it is needed before the C (i.e. if this 
            #double bond is on a side chain)
            if (is_double and 
                ((not all([(x in chainDict.keys()) for x in db_pair])) and
                 (chain_parent in db_pair)) and 
                (db_pair not in compl_dbs)):
                gecko_str += "="
                compl_dbs.append(db_pair)
            
            #add the C atom label (accounting for aromaticity)
            if root_atom.GetIsAromatic():
                gecko_str += "c"
            else:
                gecko_str += "C"
            
            #add ring label if this C is the join point for a ring
            for join_atom_i, ring_no in ringList:
                if root_atom_i == join_atom_i:
                    gecko_str += str(ring_no)
            
            #add the 'd' signifier to the double bond Cs 
            if is_double:
                gecko_str += "d"
            
            gecko_str += hDict[root_atom_i]
            
            #add the double bond symbol if it is needed after the Hs (if the 
            #double bond is in the same chain)
            if (is_double and 
                (all([(x in chainDict.keys()) for x in db_pair])) and 
                (db_pair not in compl_dbs)):
                gecko_str += "="
                compl_dbs.append(db_pair)

            for chn in chainDict[root_atom_i]:
                (ringDict2, dbList2, 
                 chainDict2, grpDict2, hDict2) = chn
                side_str, extra_atoms = GroupDict_to_GECKO(mol, ringList, 
                                                           dbList, chainDict2, 
                                                           grpDict2, hDict2,
                                                           compl_atoms, root_atom_i)
                gecko_str += f"({side_str})"
                
                for atom_i in extra_atoms:
                    compl_atoms.add(atom_i)
                
            for grp in grpDict[root_atom_i]:
                gecko_str += grp    
                
        elif root_atom.GetSymbol() == "O":
            gecko_str += "-O-"
            if hDict[root_atom_i]:
                raise Exception("Detected H atoms bound to a 'backbone' oxygen.")
            if grpDict[root_atom_i]:
                raise Exception("Detected groups bound to a 'backbone' oxygen.")
            if chainDict[root_atom_i]:
                raise Exception("Detected a side chain bound to a 'backbone' oxygen.")
        elif root_atom.GetSymbol() == "*":
            pass
        else:
            raise Exception(f"Detected a 'backbone' atom that is neither C or O: {root_atom.GetSymbol()}")
    compl_atoms.add(root_atom_i)
    return gecko_str, compl_atoms
        
def SMILES_to_GECKO(smiles):
    warnings.warn("This code is still under development and may result in incorrect GECKO strings. These will likely fail when input into GECKO, but users should double check all uses.")
    mol = SMILES_to_Mol(smiles)
    (ringDict, dbList, 
     chainDict, grpDict, hDict) = Mol_to_GroupDicts(mol)
    return GroupDict_to_GECKO(mol, ringDict, dbList, chainDict, grpDict, hDict)[0]
    

