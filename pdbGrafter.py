#!/usr/bin/env python2.7
'''
Author: Xingjie Pan
This program grafts parts of a set of structures to another set of structures.
Parameters are inputed as a sole XML file.    
'''
import sys
import collections
import Bio.PDB
import xml.etree.ElementTree as ET

class ResidueGrafter():
  '''
  residueGrafter grafts one residue from the source structure to the given position of the target structure 
  '''
  def __init__(self, a_structureT, a_modelT, a_chainT, a_residueT, a_structureS, a_modelS, a_chainS, a_residueS):
    self.structureT = a_structureT    #Target structure, type Bio.PDB.structure
    self.modelT = a_modelT            #Target model, type int
    self.chainT = a_chainT            #Target chain, type string  
    self.residueT = a_residueT        #Target residue, type int
    self.structureS = a_structureS    #Source structure, type Bio.PDB.structure
    self.modelS = a_modelS            #Source model, type int
    self.chainS = a_chainS            #Source chain, type string
    self.residueS = a_residueS        #Source residue, type int
    return
  
  def apply(self):
    '''
    Run the residue grafter
    '''
    #Get the source residue by copy
    resS = self.structureS[self.modelS][self.chainS][self.residueS].copy()
    #Change the id of the residue to the target position
    resS.id = (' ',self.residueT,' ')
    #Get the target residue
    resT = self.structureT[self.modelT][self.chainT][self.residueT]
    #Get the target residue index
    targetIndex = self.structureT[self.modelT][self.chainT].child_list.index(resT)
    #Change the residue in the list
    self.structureT[self.modelT][self.chainT].child_list[targetIndex] = resS
    #Change the residue in the dictionary
    self.structureT[self.modelT][self.chainT].child_dict[resS.id] = resS
    return


class GrafterFactory():
  '''
  grafterFactory is a wrapper of a list of grafters
  '''
  def __init__(self):
    self.grafterDic = collections.OrderedDict()
    return
  
  def add(self, a_grafterName ,a_grafter):
    '''
    Add a new grafter to the grafter factory
    '''
    if a_grafterName in self.grafterDic.keys():
      raise Exception("Grafter name "+a_grafterName+" already exists in the grafter factory!")
    else:
      self.grafterDic[a_grafterName] = a_grafter
    return 
  
  def applyAll(self):
    '''
    Run all grafters inside the grafterFactory
    '''
    for grafter in self.grafterDic.keys():
      self.grafterDic[grafter].apply()
    return


class GrafterXMLInterface():
  '''
  grafterXMLInterface handles the XML input file and create a grafterFactory to do the grafting
  '''
  def __init__(self, a_xmlName):
    self.PDBparser = Bio.PDB.PDBParser()
    self.PDBio = Bio.PDB.PDBIO()
    self.gFactory = GrafterFactory() 
    self.structureDic = {}                      #Dictionary of all structures specified in the XML file. name -> Bio.PDB.Structure
    self.outputDic = {}                         #Dictionary of structures to output. name -> file name
    self.xmlTree = ET.parse(a_xmlName)
    
    #Parsing the xml file
    xmlRoot = self.xmlTree.getroot()
    
    #Get all structures
    for structure in xmlRoot.find('structures'):
      self.structureDic[structure.find('name').text] = self.PDBparser.get_structure( structure.find('name').text, structure.find('file').text  ) 
    
    #Get all grafters
    for grafter in xmlRoot.find('grafters'): 
      if grafter.tag == 'residueGrafter':
        tST= grafter.find('target').find('structure').text
        tS = self.structureDic[ tST ]
        tMT= grafter.find('target').find('model').text
        tM = int( tMT )
        tCT= grafter.find('target').find('chain').text
        tC = tCT
        tRT= grafter.find('target').find('residue').text
        tR = int( tRT )
        sST= grafter.find('source').find('structure').text
        sS = self.structureDic[ sST ]
        sMT= grafter.find('source').find('model').text
        sM = int( sMT )
        sCT= grafter.find('source').find('chain').text
        sC = sCT
        sRT= grafter.find('source').find('residue').text 
        sR = int( sRT )
        grafterName = tST+'_'+tMT+'_'+tCT+'_'+tRT+'-'+sST+'_'+sMT+'_'+sCT+'_'+sRT
        self.gFactory.add( grafterName, ResidueGrafter(tS, tM, tC, tR, sS, sM, sC, sR) )
        
    #Get all out put structures
    for structure in xmlRoot.find('outputs'):
      self.outputDic[structure.find('name').text] = structure.find('file').text
    return
 
  def applyAll(self):
    '''
    Apply all grafters in the grafter factory
    '''
    self.gFactory.applyAll()
    return

  def output(self):
    '''
    Save pdb output files specified by the XML input file
    '''
    for structure in self.outputDic.keys():
      self.PDBio.set_structure( self.structureDic[structure] ) 
      self.PDBio.save( self.outputDic[structure] )
    return


if __name__ == '__main__':
  try:
    xmlFile = sys.argv[1] 
    gXMLI = GrafterXMLInterface(xmlFile)
    gXMLI.applyAll() 
    gXMLI.output()
  except Exception as e:
    print e
