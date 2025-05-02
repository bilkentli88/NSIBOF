import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

random_state = 1977


def create_validation_datasets(dataset_name = "Adiac",fold_number = 4):


    dataset_full_filename = "../DataSets/" + dataset_name + "_TRAIN"

    X = pd.read_csv(dataset_full_filename,sep='\s+',header=None)

    y = X[0]

    skf = StratifiedKFold(n_splits=fold_number
        , shuffle = True
        ,random_state = random_state)
    skf.get_n_splits()

    #print(skf) 


    current_fold = 1

    for train_index, test_index in skf.split(X, y):
        df_train_to_save = X.loc[train_index,:]    
        df_test_to_save = X.loc[test_index,:]    

        fold_dataset_name_train = "../DataSets/" +dataset_name + "-train-cross-validation-" + str(fold_number) + "-"  + str(current_fold) + "_TRAIN"
        fold_dataset_name_test = "../DataSets/" + dataset_name  + "-train-cross-validation-"  + str(fold_number) + "-"  + str(current_fold) + "_TEST"

        current_fold = current_fold +1
        df_train_to_save.to_csv(fold_dataset_name_train,sep=" ",header=False,index=False)
        df_test_to_save.to_csv(fold_dataset_name_test,sep=" ",header=False,index=False)

        print(fold_dataset_name_train)
        print(fold_dataset_name_test)
    #print("TRAIN:", train_index, "TEST:", test_index)

    


if __name__ == "__main__":
    import sys
    if(len(sys.argv) > 1):
        dataset_name = sys.argv[1]
        create_validation_datasets(dataset_name,4)
    else:
        create_validation_datasets("Computers",4)