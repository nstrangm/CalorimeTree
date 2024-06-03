# CalorimeTree
Analysis of isolated Photons, neutral pions, Jets and their correlations using Tree input from AliPhysics or O2Physics

# How to get started (locally):
1. Set correct file path in GJ18D/100_Charged/InputFiles/InputFiles_group_0.txt
2. Execute 'root -q -b -l Analysis/makeHistosFromTree.C+\(\"GJ18D/100_Charged/Standard\"\)'
2. Execute 'root -q -b -l Analysis/plotHistosFromTree.C+\(\"GJ18D/100_Charged/Standard\"\)'