import os
import unittest
from __main__ import vtk, qt, ctk, slicer

#
# AirwaySegmentation
#

class AirwaySegmentation:
  def __init__(self, parent):
    parent.title = "AirwaySegmentation" # TODO make this more human readable by adding spaces
    parent.categories = ["Segmentation"]
    parent.dependencies = []
    parent.contributors = ["Pietro Nardelli (University College Cork)"] 
    parent.helpText = """
    Scripted module for Airway segmentation.
    """
    parent.acknowledgementText = """
    This file was originally developed by Pietro Nardelli, University College of Cork (UCC).
""" # replace with organization, grant and thanks.
    self.parent = parent

    # Add this test to the SelfTest module's list for discovery when the module
    # is created.  Since this module may be discovered before SelfTests itself,
    # create the list if it doesn't already exist.
    try:
      slicer.selfTests
    except AttributeError:
      slicer.selfTests = {}
    slicer.selfTests['AirwaySegmentation'] = self.runTest

  def runTest(self):
    tester = AirwaySegmentationTest()
    tester.runTest()

#
# qAirwaySegmentationWidget
#

class AirwaySegmentationWidget:
  def __init__(self, parent = None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.parent.show()

  def setup(self):
    # Instantiate and connect widgets ...

    #
    # IO Area
    #
    IOCollapsibleButton = ctk.ctkCollapsibleButton()
    IOCollapsibleButton.text = "IO Parameters"
    self.layout.addWidget(IOCollapsibleButton)

    # Layout within the dummy collapsible button
    IOFormLayout = qt.QFormLayout(IOCollapsibleButton)

    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ( ("vtkMRMLScalarVolumeNode"), "" )
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = True
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    IOFormLayout.addRow("Input Volume: ", self.inputSelector)

    self.fiducialsList = slicer.qMRMLNodeComboBox()
    self.fiducialsList.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.fiducialsList.selectNodeUponCreation = False
    self.fiducialsList.addEnabled = True
    self.fiducialsList.removeEnabled = True
    self.fiducialsList.noneEnabled = False
    self.fiducialsList.showHidden = False
    self.fiducialsList.showChildNodeTypes = False
    self.fiducialsList.setMRMLScene( slicer.mrmlScene )
    self.fiducialsList.setToolTip( "Place a fiducial point within the trachea." )
    self.fiducialsList.baseName = 'AirwayFiducial'
    IOFormLayout.addRow("Seed: ", self.fiducialsList)

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Segmentation Parameters"
    parametersCollapsibleButton.setChecked(False)
    self.layout.addWidget(parametersCollapsibleButton)
    
    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)
    #
    # Label color slider
    #
    self.labelColorSliderWidget = ctk.ctkSliderWidget()
    self.labelColorSliderWidget.singleStep = 1.0
    self.labelColorSliderWidget.minimum = 1.0
    self.labelColorSliderWidget.maximum = 50.0
    self.labelColorSliderWidget.value = 2.0
    self.labelColorSliderWidget.setToolTip("Set color for the airway label")
    parametersFormLayout.addRow("Airway Label Color", self.labelColorSliderWidget)

    #
    # Apply Button
    #
    self.applyButton = qt.QPushButton("Apply")
    self.applyButton.toolTip = "Run the algorithm."
    if( self.inputSelector.currentNode() and self.fiducialsList.currentNode() ):
      self.applyButton.enabled = True
    else:
      self.applyButton.enabled = False
      
    self.applyButton.setFixedSize(150,45)
    self.layout.addWidget(self.applyButton, 0, 4)

    #
    # Link to Bronchoscopy Module
    #
    """self.bronchoscopyButton = qt.QPushButton("Link To Bronchoscopy Navigation")
    self.bronchoscopyButton.toolTip = "Connect to the Bronchoscopy module."
    #self.bronchoscopyButton.checkable = True
    self.bronchoscopyButton.enabled = False
    self.bronchoscopyButton.setFixedSize(200,50)
    self.layout.addWidget(self.bronchoscopyButton, 0, 4)"""

    # connections
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.fiducialsList.connect("currentNodeChanged(vtkMRMLNode*)",self.onSelect)
    #self.bronchoscopyButton.connect('clicked(bool)', self.onBronchoscopyButton)

    # Add vertical spacer
    self.layout.addStretch(1)

  def cleanup(self):
    pass

  def onSelect(self):
    self.applyButton.enabled = self.inputSelector.currentNode() and self.fiducialsList.currentNode()
    #self.bronchoscopyButton.enabled = False

  def onApplyButton(self):
    print("Run the algorithm")
    
    self.applyButton.enabled = False
    #self.bronchoscopyButton.enabled = False
            
    nodeType = 'vtkMRMLLabelMapVolumeNode'
    self.labelNode = slicer.mrmlScene.CreateNodeByClass(nodeType)
    self.labelNode.SetScene(slicer.mrmlScene)
    self.labelNode.SetName(slicer.mrmlScene.GetUniqueNameByString('AirwayLabel'))
    slicer.mrmlScene.AddNode(self.labelNode)

    logic = AirwaySegmentationLogic()
    labelColor = int(self.labelColorSliderWidget.value)
      
    try:
      logic.run(self.inputSelector.currentNode(), self.labelNode, self.fiducialsList.currentNode(),labelColor)
      logic.createModel(self.labelNode)
      self.applyButton.enabled = True
      #self.bronchoscopyButton.enabled = True
    except Exception, e:
      import traceback
      traceback.print_exc()
      qt.QMessageBox.warning(slicer.util.mainWindow(), 
          "Running", 'Exception!\n\n' + str(e) + "\n\nSee Python Console for Stack Trace")

  def onBronchoscopyButton(self):
    self.bronchoscopyButton.enabled = True
    mainWindow = slicer.util.mainWindow()
    mainWindow.moduleSelector().selectModule('Bronchoscopy')

#
# AirwaySegmentationLogic
#

class AirwaySegmentationLogic:
  """This class should implement all the actual 
  computation done by your module.  The interface 
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget
  """
  def __init__(self):
    pass

  def hasImageData(self,volumeNode):
    """This is a dummy logic method that 
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      print('no volume node')
      return False
    if volumeNode.GetImageData() == None:
      print('no image data')
      return False
    return True

  def delayDisplay(self,message,msec=1000):
    #
    # logic version of delay display
    #
    print(message)
    self.info = qt.QDialog()
    self.infoLayout = qt.QVBoxLayout()
    self.info.setLayout(self.infoLayout)
    self.label = qt.QLabel(message,self.info)
    self.infoLayout.addWidget(self.label)
    qt.QTimer.singleShot(msec, self.info.close)
    self.info.exec_()

  def run(self,inputVolume,outputVolume,fiducialsList,labelValue):
    """
    Run the actual algorithm
    """
    self.labelValue = labelValue

    volumeName = inputVolume.GetName()
    n = slicer.util.getNode(volumeName)
    instUIDs = n.GetAttribute('DICOM.instanceUIDs').split()
    fileName = slicer.dicomDatabase.fileForInstance(instUIDs[0])
    convolutionKernel = slicer.dicomDatabase.fileValue(fileName,'0018,1210')
    
    airwaySegmentationModule = slicer.modules.airwaysegmentationcli
    parameters = {
          "inputVolume": inputVolume.GetID(),
          "reconstructionKernelType": convolutionKernel,
          "label": outputVolume.GetID(),
          "seed": fiducialsList.GetID(),
          "labelValue": labelValue,
          }
    self.delayDisplay('Running the algorithm')
    slicer.cli.run( airwaySegmentationModule,None,parameters,wait_for_completion = True )
   
    return True
    
  def createModel(self,labelVolume):
    """
    Create the 3D model from the airway label
    """
    modelHierarchyCollection = slicer.mrmlScene.GetNodesByName('AirwayModelHierarchy')
    if( modelHierarchyCollection.GetNumberOfItems() >= 1 ):
      modelHierarchy = modelHierarchyCollection.GetItemAsObject(0)
    else:
      nodeType = 'vtkMRMLModelHierarchyNode'
      modelHierarchy = slicer.mrmlScene.CreateNodeByClass(nodeType)
      modelHierarchy.SetScene(slicer.mrmlScene)
      modelHierarchy.SetName(slicer.mrmlScene.GetUniqueNameByString('AirwayModelHierarchy'))
      slicer.mrmlScene.AddNode(modelHierarchy)
    
    parameters = {}
    parameters["InputVolume"] = labelVolume.GetID()
    parameters["ModelSceneFile"] = modelHierarchy.GetID()
    parameters["Name"] = 'AirwayModel'
    #parameters["FilterType"] = 'Laplacian'
    parameters["Smooth"] = 20
    parameters["Decimate"] = 0.10
    
    modelMaker = slicer.modules.modelmaker
    slicer.cli.run(modelMaker, None, parameters,True)
  
    lm = slicer.app.layoutManager()
    threeDView = lm.threeDWidget( 0 ).threeDView()
    threeDView.resetFocalPoint()
    threeDView.lookFromViewAxis(ctk.ctkAxesWidget().Anterior)
    
    return True

class AirwaySegmentationTest(unittest.TestCase):
  """
  This is the test case for your scripted module.
  """

  def delayDisplay(self,message,msec=1000):
    """This utility method displays a small dialog and waits.
    This does two things: 1) it lets the event loop catch up
    to the state of the test so that rendering and widget updates
    have all taken place before the test continues and 2) it
    shows the user/developer/tester the state of the test
    so that we'll know when it breaks.
    """
    print(message)
    self.info = qt.QDialog()
    self.infoLayout = qt.QVBoxLayout()
    self.info.setLayout(self.infoLayout)
    self.label = qt.QLabel(message,self.info)
    self.infoLayout.addWidget(self.label)
    qt.QTimer.singleShot(msec, self.info.close)
    self.info.exec_()

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_AirwaySegmentation1()

  def test_AirwaySegmentation1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests sould exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        print('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        print('Loading %s...\n' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading\n')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = AirwaySegmentationLogic()
    self.assertTrue( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')
