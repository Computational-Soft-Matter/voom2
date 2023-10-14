# trace generated using paraview version 5.7.0-RC1
#
# To ensure correct image size when batch processing, please search 
# for and uncomment the line `# renderView*.ViewSize = [*,*]`

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# find source
toDelete0_Psvtk = FindSource('ToDelete0_Ps.vtk')

# create a new 'Glyph'
glyph1 = Glyph(Input=toDelete0_Psvtk,
    GlyphType='Arrow')
glyph1.OrientationArray = ['POINTS', 'No orientation array']
glyph1.ScaleArray = ['POINTS', 'No scale array']
glyph1.VectorScaleMode = 'Scale by Magnitude'
glyph1.ScaleFactor = 0.16129639744758606
glyph1.GlyphTransform = 'Transform2'
glyph1.GlyphMode = 'Uniform Spatial Distribution (Bounds Based)'
glyph1.MaximumNumberOfSamplePoints = 5000
glyph1.Seed = 10339
glyph1.Stride = 1

# init the 'Arrow' selected for 'GlyphType'
glyph1.GlyphType.TipResolution = 6
glyph1.GlyphType.TipRadius = 0.1
glyph1.GlyphType.TipLength = 0.35
glyph1.GlyphType.ShaftResolution = 6
glyph1.GlyphType.ShaftRadius = 0.03
glyph1.GlyphType.Invert = 0

# init the 'Transform2' selected for 'GlyphTransform'
glyph1.GlyphTransform.Translate = [0.0, 0.0, 0.0]
glyph1.GlyphTransform.Rotate = [0.0, 0.0, 0.0]
glyph1.GlyphTransform.Scale = [1.0, 1.0, 1.0]

# find source
toDelete0vtk = FindSource('ToDelete0.vtk')

# Properties modified on glyph1
glyph1.GlyphType = 'Sphere'
glyph1.GlyphMode = 'All Points'

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')
# uncomment following to set a specific view size
# renderView1.ViewSize = [1555, 838]

# show data in view
glyph1Display = Show(glyph1, renderView1)

# trace defaults for the display properties.
glyph1Display.Representation = 'Surface'
glyph1Display.AmbientColor = [1.0, 1.0, 1.0]
glyph1Display.ColorArrayName = [None, '']
glyph1Display.DiffuseColor = [1.0, 1.0, 1.0]
glyph1Display.LookupTable = None
glyph1Display.MapScalars = 1
glyph1Display.MultiComponentsMapping = 0
glyph1Display.InterpolateScalarsBeforeMapping = 1
glyph1Display.Opacity = 1.0
glyph1Display.PointSize = 2.0
glyph1Display.LineWidth = 1.0
glyph1Display.RenderLinesAsTubes = 0
glyph1Display.RenderPointsAsSpheres = 0
glyph1Display.Interpolation = 'Gouraud'
glyph1Display.Specular = 0.0
glyph1Display.SpecularColor = [1.0, 1.0, 1.0]
glyph1Display.SpecularPower = 100.0
glyph1Display.Luminosity = 0.0
glyph1Display.Ambient = 0.0
glyph1Display.Diffuse = 1.0
glyph1Display.EdgeColor = [0.0, 0.0, 0.5]
glyph1Display.BackfaceRepresentation = 'Follow Frontface'
glyph1Display.BackfaceAmbientColor = [1.0, 1.0, 1.0]
glyph1Display.BackfaceDiffuseColor = [1.0, 1.0, 1.0]
glyph1Display.BackfaceOpacity = 1.0
glyph1Display.Position = [0.0, 0.0, 0.0]
glyph1Display.Scale = [1.0, 1.0, 1.0]
glyph1Display.Orientation = [0.0, 0.0, 0.0]
glyph1Display.Origin = [0.0, 0.0, 0.0]
glyph1Display.Pickable = 1
glyph1Display.Texture = None
glyph1Display.Triangulate = 0
glyph1Display.UseShaderReplacements = 0
glyph1Display.ShaderReplacements = ''
glyph1Display.NonlinearSubdivisionLevel = 1
glyph1Display.UseDataPartitions = 0
glyph1Display.OSPRayUseScaleArray = 0
glyph1Display.OSPRayScaleArray = 'Normals'
glyph1Display.OSPRayScaleFunction = 'PiecewiseFunction'
glyph1Display.OSPRayMaterial = 'None'
glyph1Display.Orient = 0
glyph1Display.OrientationMode = 'Direction'
glyph1Display.SelectOrientationVectors = 'None'
glyph1Display.Scaling = 0
glyph1Display.ScaleMode = 'No Data Scaling Off'
glyph1Display.ScaleFactor = 0.1774260401725769
glyph1Display.SelectScaleArray = 'None'
glyph1Display.GlyphType = 'Arrow'
glyph1Display.UseGlyphTable = 0
glyph1Display.GlyphTableIndexArray = 'None'
glyph1Display.UseCompositeGlyphTable = 0
glyph1Display.UseGlyphCullingAndLOD = 0
glyph1Display.LODValues = []
glyph1Display.ColorByLODIndex = 0
glyph1Display.GaussianRadius = 0.008871302008628845
glyph1Display.ShaderPreset = 'Sphere'
glyph1Display.CustomTriangleScale = 3
glyph1Display.CustomShader = """ // This custom shader code define a gaussian blur
 // Please take a look into vtkSMPointGaussianRepresentation.cxx
 // for other custom shader examples
 //VTK::Color::Impl
   float dist2 = dot(offsetVCVSOutput.xy,offsetVCVSOutput.xy);
   float gaussian = exp(-0.5*dist2);
   opacity = opacity*gaussian;
"""
glyph1Display.Emissive = 0
glyph1Display.ScaleByArray = 0
glyph1Display.SetScaleArray = ['POINTS', 'Normals']
glyph1Display.ScaleArrayComponent = 'X'
glyph1Display.UseScaleFunction = 1
glyph1Display.ScaleTransferFunction = 'PiecewiseFunction'
glyph1Display.OpacityByArray = 0
glyph1Display.OpacityArray = ['POINTS', 'Normals']
glyph1Display.OpacityArrayComponent = 'X'
glyph1Display.OpacityTransferFunction = 'PiecewiseFunction'
glyph1Display.DataAxesGrid = 'GridAxesRepresentation'
glyph1Display.SelectionCellLabelBold = 0
glyph1Display.SelectionCellLabelColor = [0.0, 1.0, 0.0]
glyph1Display.SelectionCellLabelFontFamily = 'Arial'
glyph1Display.SelectionCellLabelFontFile = ''
glyph1Display.SelectionCellLabelFontSize = 18
glyph1Display.SelectionCellLabelItalic = 0
glyph1Display.SelectionCellLabelJustification = 'Left'
glyph1Display.SelectionCellLabelOpacity = 1.0
glyph1Display.SelectionCellLabelShadow = 0
glyph1Display.SelectionPointLabelBold = 0
glyph1Display.SelectionPointLabelColor = [1.0, 1.0, 0.0]
glyph1Display.SelectionPointLabelFontFamily = 'Arial'
glyph1Display.SelectionPointLabelFontFile = ''
glyph1Display.SelectionPointLabelFontSize = 18
glyph1Display.SelectionPointLabelItalic = 0
glyph1Display.SelectionPointLabelJustification = 'Left'
glyph1Display.SelectionPointLabelOpacity = 1.0
glyph1Display.SelectionPointLabelShadow = 0
glyph1Display.PolarAxes = 'PolarAxesRepresentation'

# init the 'PiecewiseFunction' selected for 'OSPRayScaleFunction'
glyph1Display.OSPRayScaleFunction.Points = [0.0, 0.0, 0.5, 0.0, 1.0, 1.0, 0.5, 0.0]
glyph1Display.OSPRayScaleFunction.UseLogScale = 0

# init the 'Arrow' selected for 'GlyphType'
glyph1Display.GlyphType.TipResolution = 6
glyph1Display.GlyphType.TipRadius = 0.1
glyph1Display.GlyphType.TipLength = 0.35
glyph1Display.GlyphType.ShaftResolution = 6
glyph1Display.GlyphType.ShaftRadius = 0.03
glyph1Display.GlyphType.Invert = 0

# init the 'PiecewiseFunction' selected for 'ScaleTransferFunction'
glyph1Display.ScaleTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]
glyph1Display.ScaleTransferFunction.UseLogScale = 0

# init the 'PiecewiseFunction' selected for 'OpacityTransferFunction'
glyph1Display.OpacityTransferFunction.Points = [-0.9749279022216797, 0.0, 0.5, 0.0, 0.9749279022216797, 1.0, 0.5, 0.0]
glyph1Display.OpacityTransferFunction.UseLogScale = 0

# init the 'GridAxesRepresentation' selected for 'DataAxesGrid'
glyph1Display.DataAxesGrid.XTitle = 'X Axis'
glyph1Display.DataAxesGrid.YTitle = 'Y Axis'
glyph1Display.DataAxesGrid.ZTitle = 'Z Axis'
glyph1Display.DataAxesGrid.XTitleColor = [1.0, 1.0, 1.0]
glyph1Display.DataAxesGrid.XTitleFontFamily = 'Arial'
glyph1Display.DataAxesGrid.XTitleFontFile = ''
glyph1Display.DataAxesGrid.XTitleBold = 0
glyph1Display.DataAxesGrid.XTitleItalic = 0
glyph1Display.DataAxesGrid.XTitleFontSize = 12
glyph1Display.DataAxesGrid.XTitleShadow = 0
glyph1Display.DataAxesGrid.XTitleOpacity = 1.0
glyph1Display.DataAxesGrid.YTitleColor = [1.0, 1.0, 1.0]
glyph1Display.DataAxesGrid.YTitleFontFamily = 'Arial'
glyph1Display.DataAxesGrid.YTitleFontFile = ''
glyph1Display.DataAxesGrid.YTitleBold = 0
glyph1Display.DataAxesGrid.YTitleItalic = 0
glyph1Display.DataAxesGrid.YTitleFontSize = 12
glyph1Display.DataAxesGrid.YTitleShadow = 0
glyph1Display.DataAxesGrid.YTitleOpacity = 1.0
glyph1Display.DataAxesGrid.ZTitleColor = [1.0, 1.0, 1.0]
glyph1Display.DataAxesGrid.ZTitleFontFamily = 'Arial'
glyph1Display.DataAxesGrid.ZTitleFontFile = ''
glyph1Display.DataAxesGrid.ZTitleBold = 0
glyph1Display.DataAxesGrid.ZTitleItalic = 0
glyph1Display.DataAxesGrid.ZTitleFontSize = 12
glyph1Display.DataAxesGrid.ZTitleShadow = 0
glyph1Display.DataAxesGrid.ZTitleOpacity = 1.0
glyph1Display.DataAxesGrid.FacesToRender = 63
glyph1Display.DataAxesGrid.CullBackface = 0
glyph1Display.DataAxesGrid.CullFrontface = 1
glyph1Display.DataAxesGrid.GridColor = [1.0, 1.0, 1.0]
glyph1Display.DataAxesGrid.ShowGrid = 0
glyph1Display.DataAxesGrid.ShowEdges = 1
glyph1Display.DataAxesGrid.ShowTicks = 1
glyph1Display.DataAxesGrid.LabelUniqueEdgesOnly = 1
glyph1Display.DataAxesGrid.AxesToLabel = 63
glyph1Display.DataAxesGrid.XLabelColor = [1.0, 1.0, 1.0]
glyph1Display.DataAxesGrid.XLabelFontFamily = 'Arial'
glyph1Display.DataAxesGrid.XLabelFontFile = ''
glyph1Display.DataAxesGrid.XLabelBold = 0
glyph1Display.DataAxesGrid.XLabelItalic = 0
glyph1Display.DataAxesGrid.XLabelFontSize = 12
glyph1Display.DataAxesGrid.XLabelShadow = 0
glyph1Display.DataAxesGrid.XLabelOpacity = 1.0
glyph1Display.DataAxesGrid.YLabelColor = [1.0, 1.0, 1.0]
glyph1Display.DataAxesGrid.YLabelFontFamily = 'Arial'
glyph1Display.DataAxesGrid.YLabelFontFile = ''
glyph1Display.DataAxesGrid.YLabelBold = 0
glyph1Display.DataAxesGrid.YLabelItalic = 0
glyph1Display.DataAxesGrid.YLabelFontSize = 12
glyph1Display.DataAxesGrid.YLabelShadow = 0
glyph1Display.DataAxesGrid.YLabelOpacity = 1.0
glyph1Display.DataAxesGrid.ZLabelColor = [1.0, 1.0, 1.0]
glyph1Display.DataAxesGrid.ZLabelFontFamily = 'Arial'
glyph1Display.DataAxesGrid.ZLabelFontFile = ''
glyph1Display.DataAxesGrid.ZLabelBold = 0
glyph1Display.DataAxesGrid.ZLabelItalic = 0
glyph1Display.DataAxesGrid.ZLabelFontSize = 12
glyph1Display.DataAxesGrid.ZLabelShadow = 0
glyph1Display.DataAxesGrid.ZLabelOpacity = 1.0
glyph1Display.DataAxesGrid.XAxisNotation = 'Mixed'
glyph1Display.DataAxesGrid.XAxisPrecision = 2
glyph1Display.DataAxesGrid.XAxisUseCustomLabels = 0
glyph1Display.DataAxesGrid.XAxisLabels = []
glyph1Display.DataAxesGrid.YAxisNotation = 'Mixed'
glyph1Display.DataAxesGrid.YAxisPrecision = 2
glyph1Display.DataAxesGrid.YAxisUseCustomLabels = 0
glyph1Display.DataAxesGrid.YAxisLabels = []
glyph1Display.DataAxesGrid.ZAxisNotation = 'Mixed'
glyph1Display.DataAxesGrid.ZAxisPrecision = 2
glyph1Display.DataAxesGrid.ZAxisUseCustomLabels = 0
glyph1Display.DataAxesGrid.ZAxisLabels = []
glyph1Display.DataAxesGrid.UseCustomBounds = 0
glyph1Display.DataAxesGrid.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]

# init the 'PolarAxesRepresentation' selected for 'PolarAxes'
glyph1Display.PolarAxes.Visibility = 0
glyph1Display.PolarAxes.Translation = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.Scale = [1.0, 1.0, 1.0]
glyph1Display.PolarAxes.Orientation = [0.0, 0.0, 0.0]
glyph1Display.PolarAxes.EnableCustomBounds = [0, 0, 0]
glyph1Display.PolarAxes.CustomBounds = [0.0, 1.0, 0.0, 1.0, 0.0, 1.0]
glyph1Display.PolarAxes.EnableCustomRange = 0
glyph1Display.PolarAxes.CustomRange = [0.0, 1.0]
glyph1Display.PolarAxes.PolarAxisVisibility = 1
glyph1Display.PolarAxes.RadialAxesVisibility = 1
glyph1Display.PolarAxes.DrawRadialGridlines = 1
glyph1Display.PolarAxes.PolarArcsVisibility = 1
glyph1Display.PolarAxes.DrawPolarArcsGridlines = 1
glyph1Display.PolarAxes.NumberOfRadialAxes = 0
glyph1Display.PolarAxes.AutoSubdividePolarAxis = 1
glyph1Display.PolarAxes.NumberOfPolarAxis = 0
glyph1Display.PolarAxes.MinimumRadius = 0.0
glyph1Display.PolarAxes.MinimumAngle = 0.0
glyph1Display.PolarAxes.MaximumAngle = 90.0
glyph1Display.PolarAxes.RadialAxesOriginToPolarAxis = 1
glyph1Display.PolarAxes.Ratio = 1.0
glyph1Display.PolarAxes.PolarAxisColor = [1.0, 1.0, 1.0]
glyph1Display.PolarAxes.PolarArcsColor = [1.0, 1.0, 1.0]
glyph1Display.PolarAxes.LastRadialAxisColor = [1.0, 1.0, 1.0]
glyph1Display.PolarAxes.SecondaryPolarArcsColor = [1.0, 1.0, 1.0]
glyph1Display.PolarAxes.SecondaryRadialAxesColor = [1.0, 1.0, 1.0]
glyph1Display.PolarAxes.PolarAxisTitleVisibility = 1
glyph1Display.PolarAxes.PolarAxisTitle = 'Radial Distance'
glyph1Display.PolarAxes.PolarAxisTitleLocation = 'Bottom'
glyph1Display.PolarAxes.PolarLabelVisibility = 1
glyph1Display.PolarAxes.PolarLabelFormat = '%-#6.3g'
glyph1Display.PolarAxes.PolarLabelExponentLocation = 'Labels'
glyph1Display.PolarAxes.RadialLabelVisibility = 1
glyph1Display.PolarAxes.RadialLabelFormat = '%-#3.1f'
glyph1Display.PolarAxes.RadialLabelLocation = 'Bottom'
glyph1Display.PolarAxes.RadialUnitsVisibility = 1
glyph1Display.PolarAxes.ScreenSize = 10.0
glyph1Display.PolarAxes.PolarAxisTitleColor = [1.0, 1.0, 1.0]
glyph1Display.PolarAxes.PolarAxisTitleOpacity = 1.0
glyph1Display.PolarAxes.PolarAxisTitleFontFamily = 'Arial'
glyph1Display.PolarAxes.PolarAxisTitleFontFile = ''
glyph1Display.PolarAxes.PolarAxisTitleBold = 0
glyph1Display.PolarAxes.PolarAxisTitleItalic = 0
glyph1Display.PolarAxes.PolarAxisTitleShadow = 0
glyph1Display.PolarAxes.PolarAxisTitleFontSize = 12
glyph1Display.PolarAxes.PolarAxisLabelColor = [1.0, 1.0, 1.0]
glyph1Display.PolarAxes.PolarAxisLabelOpacity = 1.0
glyph1Display.PolarAxes.PolarAxisLabelFontFamily = 'Arial'
glyph1Display.PolarAxes.PolarAxisLabelFontFile = ''
glyph1Display.PolarAxes.PolarAxisLabelBold = 0
glyph1Display.PolarAxes.PolarAxisLabelItalic = 0
glyph1Display.PolarAxes.PolarAxisLabelShadow = 0
glyph1Display.PolarAxes.PolarAxisLabelFontSize = 12
glyph1Display.PolarAxes.LastRadialAxisTextColor = [1.0, 1.0, 1.0]
glyph1Display.PolarAxes.LastRadialAxisTextOpacity = 1.0
glyph1Display.PolarAxes.LastRadialAxisTextFontFamily = 'Arial'
glyph1Display.PolarAxes.LastRadialAxisTextFontFile = ''
glyph1Display.PolarAxes.LastRadialAxisTextBold = 0
glyph1Display.PolarAxes.LastRadialAxisTextItalic = 0
glyph1Display.PolarAxes.LastRadialAxisTextShadow = 0
glyph1Display.PolarAxes.LastRadialAxisTextFontSize = 12
glyph1Display.PolarAxes.SecondaryRadialAxesTextColor = [1.0, 1.0, 1.0]
glyph1Display.PolarAxes.SecondaryRadialAxesTextOpacity = 1.0
glyph1Display.PolarAxes.SecondaryRadialAxesTextFontFamily = 'Arial'
glyph1Display.PolarAxes.SecondaryRadialAxesTextFontFile = ''
glyph1Display.PolarAxes.SecondaryRadialAxesTextBold = 0
glyph1Display.PolarAxes.SecondaryRadialAxesTextItalic = 0
glyph1Display.PolarAxes.SecondaryRadialAxesTextShadow = 0
glyph1Display.PolarAxes.SecondaryRadialAxesTextFontSize = 12
glyph1Display.PolarAxes.EnableDistanceLOD = 1
glyph1Display.PolarAxes.DistanceLODThreshold = 0.7
glyph1Display.PolarAxes.EnableViewAngleLOD = 1
glyph1Display.PolarAxes.ViewAngleLODThreshold = 0.7
glyph1Display.PolarAxes.SmallestVisiblePolarAngle = 0.5
glyph1Display.PolarAxes.PolarTicksVisibility = 1
glyph1Display.PolarAxes.ArcTicksOriginToPolarAxis = 1
glyph1Display.PolarAxes.TickLocation = 'Both'
glyph1Display.PolarAxes.AxisTickVisibility = 1
glyph1Display.PolarAxes.AxisMinorTickVisibility = 0
glyph1Display.PolarAxes.ArcTickVisibility = 1
glyph1Display.PolarAxes.ArcMinorTickVisibility = 0
glyph1Display.PolarAxes.DeltaAngleMajor = 10.0
glyph1Display.PolarAxes.DeltaAngleMinor = 5.0
glyph1Display.PolarAxes.PolarAxisMajorTickSize = 0.0
glyph1Display.PolarAxes.PolarAxisTickRatioSize = 0.3
glyph1Display.PolarAxes.PolarAxisMajorTickThickness = 1.0
glyph1Display.PolarAxes.PolarAxisTickRatioThickness = 0.5
glyph1Display.PolarAxes.LastRadialAxisMajorTickSize = 0.0
glyph1Display.PolarAxes.LastRadialAxisTickRatioSize = 0.3
glyph1Display.PolarAxes.LastRadialAxisMajorTickThickness = 1.0
glyph1Display.PolarAxes.LastRadialAxisTickRatioThickness = 0.5
glyph1Display.PolarAxes.ArcMajorTickSize = 0.0
glyph1Display.PolarAxes.ArcTickRatioSize = 0.3
glyph1Display.PolarAxes.ArcMajorTickThickness = 1.0
glyph1Display.PolarAxes.ArcTickRatioThickness = 0.5
glyph1Display.PolarAxes.Use2DMode = 0
glyph1Display.PolarAxes.UseLogAxis = 0

# update the view to ensure updated data information
renderView1.Update()

#### saving camera placements for all active views

# current camera placement for renderView1
renderView1.CameraPosition = [4.711915337430291, -4.270561838209119, -2.1675941192149746]
renderView1.CameraFocalPoint = [0.00574648380279541, -0.01654699444770813, 0.01078149676322937]
renderView1.CameraViewUp = [0.283339779412258, 0.6665424384787116, -0.6895213898853798]
renderView1.CameraParallelScale = 1.7360177795536773

#### uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).
