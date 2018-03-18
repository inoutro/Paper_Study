import pandas as pd
import numpy as np
import csv
from scipy.stats import ttest_ind
import math
import qvalue

############################################## Data preprocessing ###############################################


# ############################################## Data Analysis ###############################################

# ################## 1. expression data 각 유전자 (row) 정규화 진행 ##################
# cols = list(expression_data_raw.columns)

# expression_data_zscore = expression_data_raw[cols[1:]]

# for i in range(1000):
# 	df2 = (expression_data_zscore.ix[i] - expression_data_zscore.ix[i].mean())/expression_data_zscore.ix[i].std(ddof=0)
# 	expression_data_zscore = expression_data_zscore.append(df2.T)

# expression_data_zscore = expression_data_zscore[1000:]
# expression_data_zscore = pd.concat([expression_data_raw[0],expression_data_zscore], axis = 1)


################## 2. pathway overlaid ##################

def pathway_overlaid(expression_data, pathway):
	df_overlaid = pd.DataFrame()
	for i in pathway:
		for j in range(len(expression_data)):
			if expression_data["GeneSymbol Name"][j] == i:
				# print(expression_data.ix[j])
				df_overlaid = df_overlaid.append(expression_data.ix[j], ignore_index=True)
				# print(df_overlaid)

	df_overlaid = df_overlaid.reset_index(drop=True)

	df_overlaid = pd.DataFrame(df_overlaid, columns = ["GeneSymbol Name", "BS16.0001439_01319724_JM", "BS16.0012807_00192199_CY", "BS16.0017904_00497033_KH", "BS16.0023523_01360069_BE", "BS16.0023862_01269828_JE", "SNUH_PA1", "SNUH_PA2", "SNUH_PA3", "SNUH_PA4", "SNUH_PA5", "SNUH_PA6", "SNUH_PA7", "SNUH_PA8", "SNUH_PA10", "SNUH_PA11"])

	return df_overlaid



################## 3. t-score / Activity apj calculate ##################

############ reordered by t-test about gene ############
def ordered_by_ttest(expression_data, label1, label2):
	sample1 = expression_data[label1]
	sample2 = expression_data[label2]

	t_value = []
	for i in range(len(expression_data)):
		ttt = ttest_ind(sample1.loc[i], sample2.loc[i])
		t_value.append(ttt[0])

	ttest_value = pd.DataFrame(t_value, columns = ["t_test"])
	expression_data_ttest = pd.concat([expression_data, ttest_value], axis = 1)

	if expression_data_ttest["t_test"].mean() >= 0:
		expression_data_ttest = expression_data_ttest.sort_values("t_test", axis = 0, ascending = True)
	else:
		expression_data_ttest = expression_data_ttest.sort_values("t_test", axis = 0, ascending = False)
	expression_data_ttest = expression_data_ttest.reset_index(drop=True)

	return expression_data_ttest


############ calculate activity score and S(Gk) ############
def activity_score_calculate(expression_data, label1, label2):
	ordered_data = ordered_by_ttest(expression_data, label1, label2)

	# calculate Activity apj
	data_len = len(ordered_data)
	ordered_data_sqrtk = ordered_data[["case", "control"]] / math.sqrt(data_len)
	compare_list = []

	# inferring S(Gk)
	for i in range(data_len):
		apj = ordered_data_sqrtk.loc[:i].sum(axis=0)
		sample1 = apj[label1]
		sample2 = apj[label2]
		sample_len = len(sample1) + len(sample2)
		score = ttest_ind(sample1, sample2)
		print(score)
		compare_list.append(score[0])
		print()
		# 첫 번째 행은 비교할 대상이 없으므로 pass
		if i == 0:
			pass
		# 두 번째 행 부터 전의 S(Gk) 값과 비교하여 S(Gk)가 더이상 증가하지 않으면 현재까지의 Activity apj 값으로 matrix 결정
		else:
			if compare_list[i-1] > compare_list[i]:
				apj_matrix = ordered_data_sqrtk.loc[:i-1].sum(axis=0)
				break
			else:
				apj_matrix = ordered_data_sqrtk.loc[:i].sum(axis=0)

	# define final matrix
	final_score = ttest_ind(apj_matrix[label1], apj_matrix[label2])
	apj_matrix["t_score"] = final_score[0]
	apj_matrix["p_value"] = final_score[1]
	apj_matrix_final = pd.DataFrame(apj_matrix)

	return apj_matrix_final

############################################## Final matrix combine each pathway activity score ###############################################

expression_data_list = ["Normalized_DESeq_0912.xlsx"]


for ex_list in expression_data_list:

	# expression data read
	expression_data_raw = pd.read_excel("./input data/" + ex_list)
	print(expression_data_raw)
	# data log
	expression_data_log = expression_data_raw[["BS16.0001439_01319724_JM", "BS16.0012807_00192199_CY", "BS16.0017904_00497033_KH", "BS16.0023523_01360069_BE", "BS16.0023862_01269828_JE", "SNUH_PA1", "SNUH_PA2", "SNUH_PA3", "SNUH_PA4", "SNUH_PA5", "SNUH_PA6", "SNUH_PA7", "SNUH_PA8", "SNUH_PA10", "SNUH_PA11"]]
	expression_data_log = np.log(expression_data_log)
	print(expression_data_log)
	expression_data_mer = pd.DataFrame(expression_data_raw, columns = ["GeneSymbol Name"])
	expression_data_merge = pd.concat([expression_data_mer, expression_data_log], axis = 1)
	print(expression_data_merge)
	print(expression_data_merge.isnull().sum())

	# pathway list read
	asd = pd.read_csv("./input data/df_kegg_pathway.txt", sep=" ")

	final_matrix = pd.DataFrame()
	pathway_list = np.unique(asd["pathway"])

	for i in pathway_list[0:1]:
		pathway_list_each = list(asd.gene[asd.pathway == i])

		expression_overlaid_pathway = pathway_overlaid(expression_data_raw, pathway_list_each)

		expression_overlaid_pathway.columns = ["gene"] + ["case"]*5 + ["control"]*10

		pathway_matrix = activity_score_calculate(expression_overlaid_pathway, "case", "control")
		pathway_matrix_T = pathway_matrix.T
		pathway_matrix_T["pathway"] = i

		final_matrix = final_matrix.append(pathway_matrix_T)

	print(final_matrix)
	print(final_matrix.loc["p_value"])
	final_matrix["q-value"] = qvalue.estimate(final_matrix["p_value"])


	print(final_matrix)

	final_matrix.to_csv("./output data/0926_final_matrix_" + ex_list.split(".")[0] + ".csv", index=False)