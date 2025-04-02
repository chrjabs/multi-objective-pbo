import gdb.printing
from prettyThing.PrettyPrintVec import StdVecPrinter
from prettyThing.PrettyPrintLit import LitPrinter
from prettyThing.PrettyPrintLbool import LboolPrinter
from prettyThing.PrettyPrintPBObjFunction import PBObjFunctionPrinter
from prettyThing.PrettyPrintPair import PairPrinter
from prettyThing.PrettyPrintModel import ModelPrinter
from prettyThing.PrettyPrintRootLits import RootLitsPrinter
from prettyThing.PrettyPrintYPoint import YPointPrinter


def build_pretty_printer():
    pp = gdb.printing.RegexpCollectionPrettyPrinter(
        "open-wbo_debug")
    pp.add_printer('vec', '^Glucose::vec<.*>$', StdVecPrinter)
    pp.add_printer('Lit', '^Glucose::Lit$', LitPrinter)
    pp.add_printer('lbool', '^Glucose::lbool$', LboolPrinter)
    pp.add_printer('pair', '^std::pair<.*,.*>$', PairPrinter)
    pp.add_printer('PBObjFunction', '^openwbo::PBObjFunction$',
                   PBObjFunctionPrinter)
    pp.add_printer('Model', '^openwbo::Model$', ModelPrinter)
    pp.add_printer('RootLits', '^rootLits::RootLits$', RootLitsPrinter)
    pp.add_printer('YPoint', '^openwbo::YPoint$', YPointPrinter)
    return pp


gdb.printing.register_pretty_printer(
    gdb.current_objfile(),
    build_pretty_printer())
