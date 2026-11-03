import cz.siret.prank.utils.PdbUtils
def structure = PdbUtils.loadFromCifFile(args[0])
def atoms = []
for (chain in structure.getModel(0)) {
    for (group in chain.getAtomGroups()) {
        for (atom in group.getAtoms()) {
            def number = group.getResidueNumber()
            atoms.add([chain.getName(), number.getSeqNum(), number.getInsCode() ?: '',
                       group.getPDBName(), atom.getName(), atom.getElement().toString(), atom.getX(), atom.getY(), atom.getZ()])
        }
    }
}
println('READ_COUNT\t'+atoms.size())
for (atom in atoms) println('ATOM\t'+atom.join('\t'))
