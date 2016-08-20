reinitialize everything;
load http://www.rcsb.org/pdb/files/1U7V.pdb;
viewport 480,480;
preset.publication("1U7V");
#show mesh;

bg_color white;

color pink, chain A; #SMAD2
color palegreen, chain B; #SMAD4
color lightblue, chain C; #SMAD2

#Residues from chain A
sele red_SMAD2_D304_A, (resi 304 and chain A); color grey, red_SMAD2_D304_A; show spheres, red_SMAD2_D304_A; #SMAD2: p.D304G, p.D274G
sele blue_SMAD2_F356_A, (resi 356 and chain A); color grey, blue_SMAD2_F356_A; show spheres, blue_SMAD2_F356_A; #SMAD2: p.F356C, p.F326C
sele blue_SMAD2_G421_A, (resi 421 and chain A); color blue, blue_SMAD2_G421_A; show spheres, blue_SMAD2_G421_A; #SMAD2: p.G421W, p.G391W
sele red_SMAD2_D450_A, (resi 450 and chain A); color red, red_SMAD2_D450_A; show spheres, red_SMAD2_D450_A; #SMAD2: p.D450N, p.D420N

#Residues from chain B
sele green_SMAD4_R361_B, (resi 361 and chain B); color grey, green_SMAD4_R361_B; show spheres, green_SMAD4_R361_B; #SMAD4: p.R361C, p.R361H, p.R361S, p.R361P, p.RFCLG361in_frame_del
sele green_SMAD4_D537_B, (resi 537 and chain B); color green, green_SMAD4_D537_B; show spheres, green_SMAD4_D537_B; #SMAD4: p.D537E, p.D537G, p.D537V, p.D537Y

#Residues from chain C
sele red_SMAD2_D304_C, (resi 304 and chain C); color grey, red_SMAD2_D304_C; show spheres, red_SMAD2_D304_C; #SMAD2: p.D304G, p.D274G
sele red_SMAD2_D450_C, (resi 450 and chain C); color red, red_SMAD2_D450_C; show spheres, red_SMAD2_D450_C; #SMAD2: p.D450N, p.D420N




#mset 1
#mdo 1: turn x,0.2; turn y,0.2; turn z,0.2;
#mplay