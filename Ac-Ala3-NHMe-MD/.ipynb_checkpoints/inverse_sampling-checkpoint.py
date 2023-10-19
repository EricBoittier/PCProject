import pandas as pd



#https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.sample.html


df = pd.read_csv("example.csv")

df.sample(n=None, frac=None, replace=False, weights=None, random_state=None, axis=None, ignore_index=False)
