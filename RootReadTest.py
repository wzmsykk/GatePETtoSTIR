import ROOT


inFile=ROOT.TFile.Open("file.root","READ")
tree=inFile.Get("Coincidences")
entries=tree.GetEntries()
print("entries:",entries)

for entry in tree:
    energy_1=tree.energy1
    print(energy_1)

