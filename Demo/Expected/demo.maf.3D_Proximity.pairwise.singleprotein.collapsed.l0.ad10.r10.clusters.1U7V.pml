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
sele green_SMAD2_S276_A, (resi 276 and chain A); color green, green_SMAD2_S276_A; show spheres, green_SMAD2_S276_A; #SMAD2: p.S276L, p.S246L
sele green_SMAD2_D304_A, (resi 304 and chain A); color green, green_SMAD2_D304_A; show spheres, green_SMAD2_D304_A; #SMAD2: p.D304G, p.D274G
sele green_SMAD2_N307_A, (resi 307 and chain A); color green, green_SMAD2_N307_A; show spheres, green_SMAD2_N307_A; #SMAD2: p.N307H, p.N277H
sele magenta_SMAD2_R321_A, (resi 321 and chain A); color grey, magenta_SMAD2_R321_A; show spheres, magenta_SMAD2_R321_A; #SMAD2: p.R321Q, p.R291Q
sele magenta_SMAD2_A323_A, (resi 323 and chain A); color magenta, magenta_SMAD2_A323_A; show spheres, magenta_SMAD2_A323_A; #SMAD2: p.A323V, p.A293V
sele brightorange_SMAD2_M327_A, (resi 327 and chain A); color brightorange, brightorange_SMAD2_M327_A; show spheres, brightorange_SMAD2_M327_A; #SMAD2: p.M327I, p.M297I
sele deepolive_SMAD2_V398_A, (resi 398 and chain A); color grey, deepolive_SMAD2_V398_A; show spheres, deepolive_SMAD2_V398_A; #SMAD2: p.V398A, p.V368A
sele green_SMAD2_L442_A, (resi 442 and chain A); color green, green_SMAD2_L442_A; show spheres, green_SMAD2_L442_A; #SMAD2: p.L442V, p.L412V
sele green_SMAD2_L446_A, (resi 446 and chain A); color green, green_SMAD2_L446_A; show spheres, green_SMAD2_L446_A; #SMAD2: p.L446V, p.L416V

#Residues from chain B
sele red_SMAD4_A327_B, (resi 327 and chain B); color red, red_SMAD4_A327_B; show spheres, red_SMAD4_A327_B; #SMAD4: p.A327V
sele chocolate_SMAD4_D355_B, (resi 355 and chain B); color chocolate, chocolate_SMAD4_D355_B; show spheres, chocolate_SMAD4_D355_B; #SMAD4: p.D355G
sele chocolate_SMAD4_P356_B, (resi 356 and chain B); color chocolate, chocolate_SMAD4_P356_B; show spheres, chocolate_SMAD4_P356_B; #SMAD4: p.P356L
sele chocolate_SMAD4_S357_B, (resi 357 and chain B); color chocolate, chocolate_SMAD4_S357_B; show spheres, chocolate_SMAD4_S357_B; #SMAD4: p.S357P
sele chocolate_SMAD4_R361_B, (resi 361 and chain B); color grey, chocolate_SMAD4_R361_B; show spheres, chocolate_SMAD4_R361_B; #SMAD4: p.R361C, p.R361S, p.R361H, p.R361P, p.RFCLG361in_frame_del
sele blue_SMAD4_A406_B, (resi 406 and chain B); color blue, blue_SMAD4_A406_B; show spheres, blue_SMAD4_A406_B; #SMAD4: p.A406T
sele blue_SMAD4_K428_B, (resi 428 and chain B); color grey, blue_SMAD4_K428_B; show spheres, blue_SMAD4_K428_B; #SMAD4: p.K428T
sele deepblue_SMAD4_R496_B, (resi 496 and chain B); color deepblue, deepblue_SMAD4_R496_B; show spheres, deepblue_SMAD4_R496_B; #SMAD4: p.R496H
sele deepblue_SMAD4_C499_B, (resi 499 and chain B); color grey, deepblue_SMAD4_C499_B; show spheres, deepblue_SMAD4_C499_B; #SMAD4: p.C499Y
sele red_SMAD4_S504_B, (resi 504 and chain B); color red, red_SMAD4_S504_B; show spheres, red_SMAD4_S504_B; #SMAD4: p.S504R
sele red_SMAD4_K519_B, (resi 519 and chain B); color red, red_SMAD4_K519_B; show spheres, red_SMAD4_K519_B; #SMAD4: p.K519N
sele red_SMAD4_W524_B, (resi 524 and chain B); color grey, red_SMAD4_W524_B; show spheres, red_SMAD4_W524_B; #SMAD4: p.W524C
sele deepblue_SMAD4_R531_B, (resi 531 and chain B); color deepblue, deepblue_SMAD4_R531_B; show spheres, deepblue_SMAD4_R531_B; #SMAD4: p.R531Q
sele deepblue_SMAD4_L533_B, (resi 533 and chain B); color deepblue, deepblue_SMAD4_L533_B; show spheres, deepblue_SMAD4_L533_B; #SMAD4: p.L533R
sele chocolate_SMAD4_D537_B, (resi 537 and chain B); color chocolate, chocolate_SMAD4_D537_B; show spheres, chocolate_SMAD4_D537_B; #SMAD4: p.D537V, p.D537G, p.D537Y, p.D537E

#Residues from chain C
sele green_SMAD2_S276_C, (resi 276 and chain C); color green, green_SMAD2_S276_C; show spheres, green_SMAD2_S276_C; #SMAD2: p.S276L, p.S246L
sele forest_SMAD2_D300_C, (resi 300 and chain C); color grey, forest_SMAD2_D300_C; show spheres, forest_SMAD2_D300_C; #SMAD2: p.D300N, p.D270N
sele magenta_SMAD2_N320_C, (resi 320 and chain C); color magenta, magenta_SMAD2_N320_C; show spheres, magenta_SMAD2_N320_C; #SMAD2: p.N320H, p.N290H
sele magenta_SMAD2_R321_C, (resi 321 and chain C); color grey, magenta_SMAD2_R321_C; show spheres, magenta_SMAD2_R321_C; #SMAD2: p.R321Q, p.R291Q
sele magenta_SMAD2_A323_C, (resi 323 and chain C); color magenta, magenta_SMAD2_A323_C; show spheres, magenta_SMAD2_A323_C; #SMAD2: p.A323V, p.A293V
sele brightorange_SMAD2_R330_C, (resi 330 and chain C); color grey, brightorange_SMAD2_R330_C; show spheres, brightorange_SMAD2_R330_C; #SMAD2: p.R330M, p.R300M
sele green_SMAD2_L442_C, (resi 442 and chain C); color green, green_SMAD2_L442_C; show spheres, green_SMAD2_L442_C; #SMAD2: p.L442V, p.L412V
sele green_SMAD2_L446_C, (resi 446 and chain C); color green, green_SMAD2_L446_C; show spheres, green_SMAD2_L446_C; #SMAD2: p.L446V, p.L416V




#mset 1
#mdo 1: turn x,0.2; turn y,0.2; turn z,0.2;
#mplay