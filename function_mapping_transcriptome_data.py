#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 17:11:54 2023

@author: sookie
"""
#%%

############### ALGORITHM HERE

#%%
def map_transcriptome_data(model,transcriptomeData):
    
    ########### NOT SUPPORTED REACTIONS
    n_noG=0
    reaction_noG = []
    for r in model.reactions:
        if str(r.gene_reaction_rule) == '':
            # print(r.id,r.gene_reaction_rule)
            n_noG+=1
            reaction_noG.append(r.id)
    
    ########### SINGLE GENE REACTIONS
    n_single=0
    reaction_single = []
    for r in model.reactions:
        if "or" not in r.gene_reaction_rule and "and" not in r.gene_reaction_rule:
            # print(r.id,r.gene_reaction_rule)
            gene_single = []
            for g in r.genes:
                gene_single.append(g)
            # print(str(gene_single))
            if len(gene_single) > 1:
                print("MISTAKE HERE !!!!!")
                break
            else:
                if str(gene_single)!="[]":
                    # print(r.id, r.gene_reaction_rule)
                    n_single+=1
                    reaction_single.append(r.id)
            
    ########### OR ONLY REACTIONS
    n_or=0
    reaction_or = []
    for r in model.reactions:
        if "or" in r.gene_reaction_rule and "and" not in r.gene_reaction_rule :
            # print(r.id, r.gene_reaction_rule)
            n_or+=1
            reaction_or.append(r.id)

    ###### AND ONLY REACTIONS
    n_and=0
    reaction_and = []

    for r in model.reactions:
        if "or" not in r.gene_reaction_rule and "and" in r.gene_reaction_rule :
            # print(r.id, r.gene_reaction_rule)
            n_and+=1
            reaction_and.append(r.id)

    ########### AND AND OR REACTIONS
    n_andOr=0
    reaction_andOr = []
    for r in model.reactions:
        if "or" in r.gene_reaction_rule and "and" in r.gene_reaction_rule :
            # print(r.id, r.gene_reaction_rule)
            n_andOr+=1
            reaction_andOr.append(r.id)
    
    ########### SUMMARY
    print("Cases no genes:",n_noG)
    print("Cases single gene:",n_single)
    print("Cases 'OR' only:",n_or)
    print("Cases 'AND' only:",n_and)
    print("Cases 'OR' and 'AND':",n_andOr)
    
    
    ########### 
    
    def length_string(gpr):
        indexes = []
        for x in range(0,len(gpr)):
            indexes.append(x)
        i = indexes[::-1]
        return i
    
    ###########
    def replace_gpr_with_fluxes(gpr,indexes,transcriptomeData):
        for x in indexes:
            # print(x,gpr[x])
            if gpr[x] == ")":
                block_indexes = []
                # print(x,gpr[x],"end here!")
                block_indexes.append(x)
                while "(" not in gpr[x]:
                    block_indexes.append(x)
                    x-=1
                    # continue
                    if gpr[x] == ")": ## another rule that has priority
                        block_indexes = []
                        # print(x,gpr[x],"end here!")
                        block_indexes.append(x)
                    # continue
                else:
                    block_indexes.append(x)
                    # print(x,gpr[x],"start here!")
                    # print("Block indexes:",block_indexes)
                    block = gpr[block_indexes[-1]:block_indexes[0]+1]
                    # print("---> Block to treat:",block)    
                    
                    ###########
                    new_flux = 0
                    ########### OR ONLY REACTIONS
                    # The OR relationship allows for alternative catalysts to each reaction; as such the total capacity is given by the sum of its components
                    if "or" in block and "and" not in block:
                        blocks = block.split("(")[1]
                        blocks = blocks.split(")")[0]
                        blocks = blocks.split(" or ")
                        # print("OR BLOCK, Genes:",blocks)
                        for g in blocks:
                            try:
                                # print(g,transcriptomeData[g])
                                new_flux = new_flux+float(transcriptomeData[g]) 
                            except:
                                new_flux = new_flux+float(g) 
                                pass
                        # print("\n--> Flux:",new_flux)
                        new_gpr = gpr.replace(block,str(new_flux))
                        # print("NEW GPR:",new_gpr)
                    ###### AND ONLY REACTIONS
                    #the maximum complex concentration is given by the minimum concentration of its components
                    if "and" in block and "or" not in block:
                        blocks = block.split("(")[1]
                        blocks = blocks.split(")")[0]
                        blocks = blocks.split(" and ")
                        # print("AND BLOCK, Genes:",blocks)
                        list_values = []
                        for g in blocks:
                            try:
                                # print(g,transcriptomeData[g])
                                list_values.append(float(transcriptomeData[g]))
                            except:
                                list_values.append(float(g))
                                pass
                        # print(list_values)
                        new_flux = min(list_values)
                        # print("\n--> Flux:",new_flux)
                        new_gpr = gpr.replace(block,str(new_flux))
                        # print("NEW GPR:",new_gpr)
                    return new_gpr      
                    break

    ###########
    reactions_to_not_contrain = []
    totToNotConst = 0
    totToNotConst_andOr = 0
    totToNotConst_and = 0
    totToNotConst_or = 0
    totToConst = 0
    reactions_rev_contraints = {}
    reactions_irrev_contraints = {}
    for r in model.reactions:
        ok=False
        new_flux = 0
        gpr = model.reactions.get_by_id(r.id).gene_reaction_rule
        old_gpr = model.reactions.get_by_id(r.id).gene_reaction_rule
        r_rev = model.reactions.get_by_id(r.id).reversibility
        # print(r.reaction)
        # print("REVERSIBILITY:",r_rev)
        all_flux_r = []
        ########### SINGLE GENE REACTIONS
        if r.id in reaction_single:
            ok=True
            for g in r.genes:
                new_flux = new_flux+float(transcriptomeData[g.id])
                all_flux_r.append(float(transcriptomeData[g.id]))
            
        
        ########### OR ONLY REACTIONS
        if r.id in reaction_or:
            ok=True 
            # The OR relationship allows for alternative catalysts to each reaction; as such the total capacity is given by the sum of its components
            for g in r.genes:
                # print(g,transcriptomeData[g.id])
                new_flux = new_flux+float(transcriptomeData[g.id])
                all_flux_r.append(float(transcriptomeData[g.id]))
                                  
        ###### AND ONLY REACTIONS
        if r.id in reaction_and:
            ok=True
            #the maximum complex concentration is given by the minimum concentration of its components
            list_values = []
            for g in r.genes:
                # print(g,transcriptomeData[g.id])
                list_values.append(float(transcriptomeData[g.id]))
                all_flux_r.append(float(transcriptomeData[g.id]))
            new_flux = min(list_values)
            
        ########### AND AND OR REACTIONS
        if r.id in reaction_andOr:
            ok=True
            gpr = model.reactions.get_by_id(r.id).gene_reaction_rule
            
            for g in r.genes:
                all_flux_r.append(float(transcriptomeData[g.id]))
                
            while ")" in str(gpr) or "(" in str(gpr):
                # print("WHILE")
                indexes = length_string(gpr)
                gpr = replace_gpr_with_fluxes(gpr,indexes,transcriptomeData)
            else:
                ########### OR ONLY REACTIONS LEFT
                # The OR relationship allows for alternative catalysts to each reaction; as such the total capacity is given by the sum of its components
                if "or" in str(gpr) and "and" not in str(gpr):
                    f=gpr.split(" or ")
                    for g in f:
                        try:
                            new_flux = new_flux+float(transcriptomeData[g])
                        except:
                            new_flux = new_flux+float(g)
                ###### AND ONLY REACTIONS LEFT
                #the maximum complex concentration is given by the minimum concentration of its components
                elif "and" in str(gpr) and "or" not in str(gpr):
                    f=gpr.split(" and ")
                    list_values = []
                    for g in f:
                        try:
                            list_values.append(float(transcriptomeData[g]))
                        except:
                            list_values.append(float(g))
                    new_flux = min(list_values)
                else:
                    print("-------------------\n")
                    print("-------------------\n")
                    print("--------------------\n")
                    print(r.id,model.reactions.get_by_id(r.id).gene_reaction_rule)
                    print(gpr)
                    print("CRAZY CASE")
                    break

        for val in all_flux_r:
            if float(val) < 10:  
                # print(r.id,"TO NOT CONSTRAIN")
                ok = False
        ########### RESULTS
        
        if ok :
            # print("\n---",r.id,"GPR:",old_gpr)
            # print("--GPR:",gpr)
            # print("--Flux:",new_flux)
            if r_rev == True:
                # print(r.reaction)
                reactions_rev_contraints[r.id] = float(new_flux)
                r.lower_bound = -float(new_flux)
                r.upper_bound = float(new_flux)
            else:
                # print(r.reaction)
                # print("NOO")
                reactions_irrev_contraints[r.id] = float(new_flux)
                r.upper_bound = float(new_flux)
            totToConst+=1
        else:
            totToNotConst+=1
            if r.id in reaction_andOr:
                totToNotConst_andOr+=1
            if r.id in reaction_or:
                totToNotConst_or+=1
            if r.id in reaction_and:
                totToNotConst_and+=1
    
    all_constraints = {**reactions_rev_contraints, **reactions_irrev_contraints}
    max_value = max(all_constraints.values())
    for k,v in all_constraints.items():
        new_v=v*1000/max_value
        all_constraints[k]=new_v
    
    for k,v in all_constraints.items():
        if k in reactions_rev_contraints.keys():
            r=model.reactions.get_by_id(k)
            r.lower_bound = -float(v)
            r.upper_bound = float(v)
        elif k in reactions_irrev_contraints.keys():
            r=model.reactions.get_by_id(k)
            r.upper_bound = float(v)

    print("Total number of reactions that will be constrained:",totToConst)
    print("Total number of reactions that won't be constrained:",totToNotConst,"(And:",totToNotConst_and,"; Or:",totToNotConst_or,"; AndOr:",totToNotConst_andOr,")")

    return model
#%%




###### MANUAL VERIFICATION COMPLEX CASE WITH "AND" + "OR" in GPR
#transcriptomeData = dico_22degrees_wort
#r="r_0001"
#gpr = model.reactions.get_by_id(r).gene_reaction_rule    
#
#old_gpr = gpr
#y = 1
#while ")" in str(gpr) or "(" in str(gpr):
#    # print("WHILE")
#    indexes = length_string(gpr)
#    gpr = replace_gpr_with_fluxes(gpr,indexes,transcriptomeData)  
#    # print("\n---Round",y,"NEW:\n",gpr) 
#    # print("\n-------OLD:\n",old_gpr)
#    # y+=1
## else:
##     print("NEW GPR:",gpr)
#print(gpr)
#
### (( 483.071718559283 + 41.2554017387574) and (1.93349209085493 + 1.69429647902912)) 
### or ((483.071718559283 + 41.2554017387574) and (235.913178727607 + 411.683270343694)) 
### or ((1.93349209085493 + 1.69429647902912) and (5217.96742810375 + 284.127543338199)) 
### or ((5217.96742810375 + 284.127543338199) and (235.913178727607 + 411.683270343694))
#
### (524.327120298 and 3.62778857) 
### or (524.327120298 and 647.596449071) 
### or (3.62778857 and 5502.094971442) 
### or (5502.094971442 and 647.596449071)
#
### (524.327120298 and [3.62778857]) 
### or ([524.327120298] and 647.596449071) 
### or ([3.62778857] and 5502.094971442) 
### or (5502.094971442 and [647.596449071])
#
### 3.62778857 + 524.327120298 + 3.62778857 + 647.596449071
### TOTAL: 1179.179146509
#
#print(new_model.reactions.get_by_id(r).lower_bound)
#print(new_model.reactions.get_by_id(r).upper_bound)
#
##%%
#print(reaction_andOr)
##%%
#transcriptomeData=dico_22degrees_SD
#
#for r in reaction_andOr:
#    ok = True
#    for g in model.reactions.get_by_id(r).genes:
#        if float(transcriptomeData[g.id]) < 10:
#            ok = False
#    if ok:
#        print("\n-",r,model.reactions.get_by_id(r).name)
#        print(model.reactions.get_by_id(r).gene_reaction_rule)
#        for g in model.reactions.get_by_id(r).genes:
#            print(g.id,g.name,transcriptomeData[g.id])
##%%
###### MANUAL VERIFICATION COMPLEX CASE WITH "AND" + "OR" in GPR
#transcriptomeData = dico_22degrees_wort
#r="r_0439"
#gpr = model.reactions.get_by_id(r).gene_reaction_rule    
#
#old_gpr = gpr
#print("GPR:",old_gpr)
#y = 1
#while ")" in str(gpr) or "(" in str(gpr):
#    # print("WHILE")
#    indexes = length_string(gpr)
#    gpr = replace_gpr_with_fluxes(gpr,indexes,transcriptomeData)  
#    # print("\n---Round",y,"NEW:\n",gpr) 
#    # print("\n-------OLD:\n",old_gpr)
#    # y+=1
## else:
##     print("NEW GPR:",gpr)
#print(gpr)
#
#for g in model.reactions.get_by_id(r).genes:
#    print(g.id,transcriptomeData[g.id])
#
### ((SPGP0DC00630 or SPGP0DN00680) and SPGP0Y01770 and (SPGP0I00570 or SPGP0E00490) and (SPGP0I00430 or SPGP0E00370) 
### and SPGP0AH00350 and SPGP0F01050 and SPGP0V00530 and (SPGP0G02940 or SPGP0R00500) and (SPGP0K02220 or SPGP0L00690) and SPGP0AC01080) 
### or 
### ((SPGP0DC00630 or SPGP0DN00680) and SPGP0Y01770 and (SPGP0I00570 or SPGP0E00490) and SPGP0AH00350 and SPGP0F01050 
### and SPGP0V00530 and (SPGP0G02940 or SPGP0R00500) and SPGP0G01100 and (SPGP0K02220 or SPGP0L00690) and SPGP0AC01080)
#
#
### ((152.424282657451 or 45.8906343673884) and 269.80106979075 and (3.47011619385828 or 5.96273815114047) and (1.93349209085493 or 1.69429647902912) 
### and 371.405431022875 and 115.626454947252 and 95.5179705240057 and (0.483373022713733 or 5177.40515836488) and (1209.65850674892 or 220.654939114733) and 64.2230930162529) 
### or 
### ((152.424282657451 or 45.8906343673884) and 269.80106979075 and (3.47011619385828 or 5.96273815114047) and 371.405431022875 and 115.626454947252 
### and 95.5179705240057 and (0.483373022713733 or 5177.40515836488) and 235.913178727607 and (1209.65850674892 or 220.654939114733) and 64.2230930162529)
#
#
### ((198.314917025) and 269.80106979075 and (9.432854345) and (3.62778857) 
### and 371.405431022875 and 115.626454947252 and 95.5179705240057 and (5177.888531388) and (1430.313445864) and 64.2230930162529)
#    
### or 
#    
### ((198.314917025) and 269.80106979075 and (9.432854345) and 371.405431022875 and 115.626454947252 
### and 95.5179705240057 and (5177.888531388) and 235.913178727607 and (1430.313445864) and 64.2230930162529)
#
### 3.62778857 + 9.432854345
#    
### TOTAL: 13.060642915
#
#print(new_model.reactions.get_by_id(r).lower_bound)
#print(new_model.reactions.get_by_id(r).upper_bound)
#
##%%
#print(reaction_and)
##%%
###### MANUAL VERIFICATION CASE WITH "AND"
#transcriptomeData = dico_22degrees_wort
#r="r_3507"
#gpr = model.reactions.get_by_id(r).gene_reaction_rule    
#print(r,gpr)
#
### SPGP0Y00260 and SPGP0CZ02130
#print("SPGP0Y00260",transcriptomeData["SPGP0Y00260"])
#print("SPGP0CZ02130",transcriptomeData["SPGP0CZ02130"])
#
### TOTAL: 76.5057299456432
#
#print(new_model.reactions.get_by_id(r).lower_bound)
#print(new_model.reactions.get_by_id(r).upper_bound)
#
##%%
#print(reaction_or)
##%%
###### MANUAL VERIFICATION CASE WITH "AND"
#transcriptomeData = dico_22degrees_wort
#r="r_4705"
#gpr = model.reactions.get_by_id(r).gene_reaction_rule    
#print(r,gpr)
#
### SPGP0DD00140 or SPGP0B00110
#print("SPGP0DD00140",transcriptomeData["SPGP0DD00140"])
#print("SPGP0B00110",transcriptomeData["SPGP0B00110"])
### 42.0067433059312 + 253.737688627161
### TOTAL: 295.744431933
#
#print(new_model.reactions.get_by_id(r).lower_bound)
#print(new_model.reactions.get_by_id(r).upper_bound)
