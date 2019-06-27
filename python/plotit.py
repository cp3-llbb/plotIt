"""
plotIt using matplotlib

Based on https://github.com/cp3-llbb/plotIt

WARNING: very much work-in-progress, many things are not implemented yet
"""

from itertools import chain
from collections import OrderedDict as odict

class BaseYAMLObject(object):
    """ Base class for objects constructed from YAML """
    required_attributes = set()
    optional_attributes = dict()
    @staticmethod
    def _normalizeAttr(att):
        return att.replace("-", "_")
    @staticmethod
    def _normalizeKeys(aDict):
        return dict((BaseYAMLObject._normalizeAttr(k), v) for k, v in aDict.iteritems())
    def __init__(self, **kwargs):
        if not all( nm in kwargs for nm in self.__class__.required_attributes ):
            raise KeyError("Attribute(s) {1} required for class {0} but not found (full dictionary: {2})".format(self.__class__.__name__, ", ".join("'{}'".format(k) for k in self.__class__.required_attributes if k not in kwargs), str(kwargs)))
        if not all( k in self.__class__.required_attributes or k in self.__class__.optional_attributes for k in kwargs.iterkeys() ):
            raise KeyError("Unknown attribute(s) for class {0}: {1} (all attributes: {2})".format(self.__class__.__name__, ", ".join("'{}'".format(k) for k in kwargs.iterkeys() if k not in self.__class__.required_attributes and k not in self.__class__.optional_attributes), str(kwargs)))
        self.__dict__.update(BaseYAMLObject._normalizeKeys(self.__class__.optional_attributes))
        self.__dict__.update(BaseYAMLObject._normalizeKeys(kwargs))
    def __repr__(self):
        return "{0}({1})".format(self.__class__.__name__, ", ".join("{0}={1!r}".format(k,getattr(self, BaseYAMLObject._normalizeAttr(k))) for k in chain(self.__class__.required_attributes, (dk for dk,dv in self.__class__.optional_attributes.iteritems() if dv is not getattr(self, BaseYAMLObject._normalizeAttr(dk))))))

def mergeDicts(first, second):
    """ trivial helper: copy first, update with second and return """
    md = dict(first)
    md.update(second)
    return md

class PlotStyle(BaseYAMLObject):
    optional_attributes = {
              "legend"          : ""
            , "legend-style"    : None
            , "legend-order"    : 0
            , "drawing-options" : ""
            , "marker-size"     : None
            , "marker-color"    : None
            , "marker-type"     : None
            , "fill-color"      : None
            , "fill-type"       : None
            , "line-width"      : None
            , "line-color"      : None
            , "line-type"       : None
            }
    def __init__(self, **kwargs):
        typ = kwargs.get("type", "MC").upper()
        super(PlotStyle, self).__init__(**kwargs)

        if typ == "MC":
            self.legend_style = "lf"
        elif typ == "SIGNAL":
            self.legend_style = "l"
        elif typ == "MC":
            self.legend_style = "pe"

        if not self.drawing_options:
            if typ in ("MC", "SIGNAL"):
                self.drawing_options = "hist"
            elif typ == "DATA":
                self.drawing_options = "pe"

        if self.fill_color:
            self.fill_color = PlotStyle.loadColor(self.fill_color)
        if self.line_color:
            self.line_color = PlotStyle.loadColor(self.line_color)
        if self.marker_color:
            self.marker_color = PlotStyle.loadColor(self.marker_color)

    @staticmethod
    def loadColor(color):
        if not ( color[0] == "#" and ( len(color) == 7 or len(color) == 9 ) ):
            raise ValueError("Color should be in the format '#RRGGBB' or '#RRGGBBAA' (got {!r})".format(color))
        r = int(color[1:3], base=16)/255.
        g = int(color[3:5], base=16)/255.
        b = int(color[5:7], base=16)/255.
        a = int(color[7:9], base=16)/255. if len(color) > 7 else 1.
        return (r,g,b,a)


from .systematics import HistoKey

class HistoFile(PlotStyle):
    required_attributes = set(("path", "type"))
    optional_attributes = mergeDicts(PlotStyle.optional_attributes, {
              "pretty-name"      : None
            , "type"             : "MC"
            ##
            , "scale"            : 1.
            , "cross-section"    : None
            , "branching-ratio"  : 1.
            , "generated-events" : None
            ##
            , "group"            : None ## ??
            , "yields-title"     : None ## ??
            , "yields-group"     : None
            , "legend-group"     : None
            ##
            , "order"            : None
            })
    def __init__(self, **kwargs):
        super(HistoFile, self).__init__(**kwargs)
        if self.pretty_name is None:
            self.pretty_name = self.path
        if self.yields_group is None:
            if self.group is not None:
                self.yields_group = self.group
            elif self.legend is not None:
                self.yields_group = self.legend
            else:
                self.yields_group = self.path ## FIXME path

        #import os.path
        #if not os.path.isfile(self.path):
        #    raise RuntimeError("File '{}' does not exist (or is not a regular file)".format(self.path))
        from cppyy import gbl
        self._tf = gbl.TFile.Open(self.path, "READ")
        #if not self._tf:
        #    raise RuntimeError("Could not open file '{}'".format(self.path))
    def getKey(self, name, **kwargs):
        return HistoKey(self._tf, name, **kwargs)

class Point(BaseYAMLObject):
    required_attributes = set(("x", "y"))
    def __init__(self, **kwargs):
        super(Point, self).__init__(**kwargs)

class Label(BaseYAMLObject):
    required_attributes = set(("text", "position"))
    optional_attributes = { "size" : 18 }
    def __init__(self, **kwargs):
        super(Label, self).__init__(**kwargs)
        self.position = Point(self.position)

class Plot(BaseYAMLObject):
    required_attributes = set(("name",))
    optional_attributes = {
              "exclude"                   : ""
            , "x-axis"                    : ""
            , "x-axis-format"             : ""
            , "y-axis"                    : ""
            , "y-axis-format"             : None
            , "normalized"                : False
            , "no-data"                   : False
            , "override"                  : False
            , "log-y"                     : False
            , "log-x"                     : False
            , "save-extensions"           : tuple()
            , "book-keeping-folder"       : ""
            , "show-ratio"                : False
            , "fit"                       : False
            , "fit-ratio"                 : False
            #
            , "fit-function"              : ""
            , "fit-legend"                : ""
            , "fit-legend-position"       : None
            , "fit-range"                 : None
            #
            , "ratio-fit-function"        : ""
            , "ratio-fit-legend"          : ""
            , "ratio-fit-legend-position" : None
            , "ratio-fit-range"           : None
            #
            , "show-errors"               : True
            , "x-axis-range"              : None
            , "y-axis-range"              : None
            , "ratio-y-axis-range"        : None
            , "blinded-range"             : None
            , "y-axis-show-zero"          : None
            , "inherits-from"             : None
            , "rebin"                     : 1
            , "labels"                    : []
            , "extra-label"               : None
            , "legend-position"           : None
            , "legend-columns"            : None
            , "show-overflow"             : None
            , "errors-type"               : None
            #
            , "binning-x"                 : None
            , "binning-y"                 : None
            , "draw-string"               : None
            , "selection-string"          : None
            #
            , "for-yields"                : True
            , "yields-title"              : True
            , "yields-table-order"        : True
            , "sort-by-yields"            : False
            #
            , "vertical-lines"            : []
            , "horizontal-lines"          : []
            , "lines"                     : []
            }
    def __init__(self, **kwargs):
        super(Plot, self).__init__(**kwargs)
        #if self.x_axis_range is not None:
        #    try:
        #        lims = tuple(float(tok.strip()) for tok in self.x_axis_range.strip("[]").split(","))
        #        if len(lims) != 2:
        #            raise ValueError("")
        #    except Exception, e:
        #        raise ValueError("Could not parse x-axis-range {0}: {1}".format(self.x_axis_range, e))
        #    self.x_axis_range = lims
        self.labels = [ Label(lblNd) for lblNd in self.labels ]

def _plotit_loadWrapper(fpath):
    """ yaml.safe_load from path """
    import yaml
    with open(fpath) as f:
        res = yaml.safe_load(f)
    return res

import os.path

def _load_includes(cfgDict, basePath):
    updDict = dict()
    for k,v in cfgDict.iteritems():
        if isinstance(v, dict):
            if len(v) == 1 and next(v.iterkeys()) == "include":
                vals = v[next(v.iterkeys())]
                newDict = dict()
                for iPath in vals:
                    if not os.path.isabs(iPath):
                        iPath = os.path.join(basePath, iPath)
                    if not os.path.exists(iPath):
                        raise IOError("Included path '{}' does not exist".format(iPath))
                    newDict.update(_plotit_loadWrapper(iPath))
                updDict[k] = newDict
                _load_includes(newDict, basePath)
            else:
                _load_includes(v, basePath)
    cfgDict.update(updDict)

def makeSystematic(item):
    from .systematics import ShapeSystVar, ConstantSystVar, LogNormalSystVar
    if isinstance(item, str):
        return ShapeSystVar(item)
    elif isinstance(item, dict):
        if len(item) == 1:
            name, val = next(item.iteritems())
            if isinstance(val, float):
                return ConstantSystVar(name, val)
            elif isinstance(val, dict):
                if val["type"] == "shape":
                    syst = ShapeSystVar(name)
                elif val["type"] == "const":
                    syst = ConstantSystVar(name, val["value"])
                elif val["type"] in ("lognormal", "ln"):
                    cfg = dict((k.replace("-", "_"), v) for k, v in val.iteritems())
                    cfg.pop("type")
                    prior = cfg.pop("prior")
                    syst = LogNormalSystVar(name, prior, **cfg)
                if "pretty-name" in val:
                    syst.pretty_name = val["pretty-name"]
                if True in val:
                    import re
                    pat = re.compile(val[True])
                    syst.on = ( lambda aPat : ( lambda fName,fObj : bool(aPat.match(fName)) ) )(pat) ## FIXME on is automatically parsed to True
                return syst
    else:
        raise ValueError("Invalid systematics node, must be either a string or a map")

def _plotIt_histoPath(histoPath, cfgRoot, baseDir):
    if os.path.isabs(histoPath):
        return histoPath
    elif os.path.isabs(cfgRoot):
        return os.path.join(cfgRoot, histoPath)
    else:
        return os.path.join(baseDir, cfgRoot, histoPath)

def plotIt_load(mainPath, histoBaseDir):
    ## load config, with includes
    cfg = _plotit_loadWrapper(mainPath)
    basedir = os.path.dirname(mainPath)
    _load_includes(cfg, basedir)
    plotDefaults = dict((k,v) for k,v in cfg["configuration"].iteritems() if k in ("y-axis-format", "show-overflow", "errors-type"))
    ## construct objects
    files = odict(sorted(dict((k, HistoFile(path=_plotIt_histoPath(k, cfg["configuration"]["root"], histoBaseDir), **v)) for k, v in cfg["files"].iteritems()).iteritems(), key=lambda (k,v) : v.order))
    ## TODO groups
    plots = dict((k, Plot(name=k, **mergeDicts(plotDefaults, v))) for k, v in cfg.get("plots", {}).iteritems())
    systematics = [ makeSystematic(item) for item in cfg.get("systematics", []) ]

    return cfg, files, plots, systematics

def getScaleForFile(f, config):
    """ Infer the scale factor for histograms from the file dict and the overall config """
    if f.type == "data":
        return 1.
    else:
        mcScale = ( config["luminosity"]*f.cross_section*f.branching_ratio / f.generated_events )
        if config.get("ignore-scales", False):
            return mcScale
        else:
            return mcScale*config.get("scale", 1.)*f.scale

def drawPlot(plot, expStack, obsStack, outdir="."):
    from .histstacksandratioplot import THistogramRatioPlot
    theplot = THistogramRatioPlot(expected=expStack, observed=obsStack) ## TODO more opts?
    theplot.draw()
    #
    if plot.x_axis_range:
        theplot.ax.set_xlim(*plot.x_axis_range)
    if plot.x_axis:
        theplot.rax.set_xlabel(plot.x_axis)
    #
    if plot.y_axis_range:
        theplot.ax.set_ylim(*plot.y_axis_range)
    else:
        if not plot.log_y:
            theplot.ax.set_ylim(0.)
    if plot.y_axis:
        theplot.ax.set_ylabel(plot.y_axis)
    elif plot.y_axis_format:
        pass
    #
    import os.path
    for ext in plot.save_extensions:
        theplot.fig.savefig(os.path.join(outdir, "{0}.{1}".format(plot.name, ext)))

def plotIt(plots, files, systematics=None, config=None):
    ## default kwargs
    if systematics is None:
        systematics = list()
    if config is None:
        config = dict()

    scaleAndSystematicsPerFile = odict((f,
        (getScaleForFile(f, config), dict((syst.name, syst) for syst in systematics if syst.on(fN, f)))
        ) for fN,f in files.iteritems())

    from .histstacksandratioplot import THistogramStack
    from .systematics import SystVarsForHist
    for pName, aPlot in plots.iteritems():
        obsStack = THistogramStack()
        expStack = THistogramStack()
        for f, (fScale, fSysts) in scaleAndSystematicsPerFile.iteritems():
            hk = f.getKey(pName, scale=fScale, rebin=aPlot.rebin, xOverflowRange=(aPlot.x_axis_range if aPlot.show_overflow else None))
            if f.type == "data":
                obsStack.add(hk, systVars=SystVarsForHist(hk, fSysts)) ##, label=..., drawOpts=...
            elif f.type == "mc":
                expStack.add(hk, systVars=SystVarsForHist(hk, fSysts), drawOpts={"fill_color":f.fill_color}) ##, label=..., drawOpts=...

        drawPlot(aPlot, expStack, obsStack)


def plotItFromYAML(yamlFileName, histoBaseDir):
    cfg, files, plots, systematics = plotIt_load(yamlFileName, histoBaseDir)
    ### get list of files, get list of systs, dict of systs per file; then list of plots: for each plot build the stacks and draw
    ## TODO cfg -> config
    plotIt(plots, files, systematics=systematics, config=cfg["configuration"])

if __name__ == "__main__": ## quick test of basic functionality
    import ROOT
    ROOT.PyConfig.IgnoreCommandLineOptions = True
    import os.path
    my_plotit_dir = "" ## FIXME
    plotItFromYAML(os.path.join(my_plotit_dir, "examples/example.yml"))
    from matplotlib import pyplot as plt
    plt.show()
