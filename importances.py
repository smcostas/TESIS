# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 13:10:51 2024

@author: santi
"""
import os 
os.chdir('OneDrive/Escritorio/PC-IMBIV/doctorado/asteraceae - muestreo/segundo capitulo/tablas segundo capitulo')
from sklearn.tree import DecisionTreeRegressor, plot_tree
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.preprocessing import OneHotEncoder
#import seaborn as sns
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.ensemble import RandomForestRegressor
from sklearn.inspection import permutation_importance
#%%
df = pd.read_csv('df_ev.csv')
df[df['Tribu'] == 'Coreopsideae']

#%%

s_model = DecisionTreeRegressor(random_state=0) ### voy a dejarlo por defecto
#%%
### pipelinde trabajo
y = df['sv']
X = df.iloc[:,1:]

## one hot encoder
types = X.dtypes
cols_cat  = types[types == 'object'].index.to_list()
cols_num  = types[types != 'object'].index.to_list()
X[cols_cat].info() ## 

encoder = OneHotEncoder(handle_unknown='ignore', sparse_output=False)
X_cat = encoder.fit_transform(X[cols_cat])
new_columns = []
for col, col_values in zip(cols_cat, encoder.categories_):
  for col_value in col_values:
    new_columns.append('{}={}'.format(col, col_value))

X = np.hstack([X_cat, X[cols_num].values])
new_columns.extend(cols_num)

X =  pd.DataFrame(data=X, columns=new_columns)
#%%%
##fiteo el modeloo
s_model.fit(X,y) ## uso completo por que quiero hacer seleccion de variables
#%% Ploteo del tree
fig = plt.figure(figsize=(25,20))
 #add title to the plot
fig.suptitle(' Tree Regressor\n selecting variables', fontsize=50)
_ = plot_tree(s_model,
                    feature_names=X.columns,
                    filled=True)


fig = plt.figure(figsize=(50, 30))
_ = plot_tree(
    s_model,
    feature_names=X.columns,
    filled=True,
    fontsize=15,  # Cambia el tamaño del texto según lo necesites
    max_depth = 5# Hace que los cuadros sean redondeados para mejorar la legibilidad
)


plt.show()


#%% importancia de cada variable
importances = s_model.feature_importances_

feature_importances = pd.DataFrame({
    'Feature': X.columns,
    'Importance': importances
}).sort_values(by='Importance', ascending=False)


## tengo que reunir las variables dummies en  originales
## primero mapeamos las variables dummies a sus originales
feature_importances['Original Feature'] = feature_importances['Feature'].apply(lambda x: x.split('=')[0] if '=' in x else x) ## spliteo el nombre de la columna y me quedo con el primer termino (el nombre original), cuando no incluye '=' es una variable numerica, se copia el mismo nombre
### ahora agrupamos
grouped_importances = feature_importances.groupby('Original Feature').sum().reset_index()
## ordeno por importancia
grouped_importances = grouped_importances.sort_values(by='Importance', ascending=False)


#%% ploteo
plt.figure(figsize=(18, 18))
plt.barh(grouped_importances['Original Feature'], grouped_importances['Importance'], color='skyblue')
plt.xlabel('Importance', fontsize = 20)
plt.xticks(fontsize=15) 
plt.yticks(fontsize = 20)
plt.gca().invert_yaxis()  # To display the most important features at the top
plt.grid(visible = True, axis = 'x') 
plt.show()

## ploteando las dummies
plt.figure(figsize=(18, 18))
plt.barh(feature_importances['Feature'], feature_importances['Importance'], color='skyblue')
plt.xlabel('Importance', fontsize = 20)
plt.gca().invert_yaxis()  # To display the most important features at the top
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)
plt.grid(visible = True, axis = 'x') 
plt.show()
#%%% metricas

preds = s_model.predict(X)

mse = mean_squared_error(y, preds)

r2 = r2_score(y, preds)
print(f'MSE: {mse}')  ## 4.0665872296158213e-17
print(f'R-squared: {r2}') ## 0.9999295718481932
## esta sobre ajustado que es lo que esperamos. para poder mostrar la importancia de las variables

#%%% USANDO RANDOM FOREST

model = RandomForestRegressor(n_estimators=100, random_state=42)
X = X.drop('pl_m', axis = 1)
model.fit(X, y)

perm_importance = permutation_importance(model, X, y, n_repeats=10, random_state=42, scoring='neg_mean_squared_error')

# Calcular la importancia total
total_importance = perm_importance.importances_mean.sum()

# Calcular la importancia relativa en porcentaje
feature_importances = pd.DataFrame({
    'Feature': X.columns,
    'Importance': perm_importance.importances_mean,
    'Importance (%)': (perm_importance.importances_mean / total_importance) * 100
}).sort_values(by='Importance', ascending=False)

feature_importances['Original Feature'] = feature_importances['Feature'].apply(lambda x: x.split('=')[0] if '=' in x else x)


grouped_importances = feature_importances.groupby('Original Feature').sum().reset_index()
grouped_importances = grouped_importances.sort_values(by='Importance (%)', ascending=False)
# Mostrar los resultados
print(feature_importances)
#%%%
# Visualización de la importancia en porcentaje
feature_importances = feature_importances[feature_importances['Importance (%)'] > 0.21]
plt.figure(figsize=(18, 18))
plt.barh(feature_importances['Feature'], feature_importances['Importance (%)'], color='skyblue')
plt.xlabel('Importancia de permutación (%)', fontsize = 20)
plt.xticks([0,20,40,60],fontsize=15) 
plt.xlim(left=0)
plt.yticks(fontsize = 20)
plt.gca().invert_yaxis()
plt.grid(visible = True, axis = 'x')
plt.show()

# visualizacion de las variables originales
plt.figure(figsize=(18, 18))
plt.barh(grouped_importances['Original Feature'], grouped_importances['Importance (%)'], color='skyblue')
plt.xlabel('Importance (%)', fontsize = 20)
plt.xticks(fontsize=15) 
plt.yticks(fontsize = 20)
plt.gca().invert_yaxis()  # To display the most important features at the top
plt.grid(visible = True, axis = 'x') 
plt.show()

#%%%
### 
sub_df = df.copy()
sub_df = sub_df[sub_df['tipo_papus'] == 'setose']
sub_df = sub_df[sub_df['Tribu'] != 'Barnadesieae']
sub_y = sub_df['sv']
sub_X = sub_df.drop(['tipo_papus','sv'], axis = 1 )


types = sub_X.dtypes
cols_cat  = types[types == 'object'].index.to_list()
cols_num  = types[types != 'object'].index.to_list()

encoder = OneHotEncoder(handle_unknown='ignore', sparse_output=False)
X_cat = encoder.fit_transform(sub_X[cols_cat])
new_columns = []
for col, col_values in zip(cols_cat, encoder.categories_):
  for col_value in col_values:
    new_columns.append('{}={}'.format(col, col_value))

sub_X = np.hstack([X_cat, sub_X[cols_num].values])
new_columns.extend(cols_num)

sub_X =  pd.DataFrame(data=sub_X, columns=new_columns)


#%% random forest setosas


model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(sub_X, sub_y)

perm_importance = permutation_importance(model, sub_X, sub_y, n_repeats=100, random_state=42, scoring='neg_mean_squared_error')

# Calcular la importancia total
total_importance = perm_importance.importances_mean.sum()

# Calcular la importancia relativa en porcentaje
feature_importances = pd.DataFrame({
    'Feature': sub_X.columns,
    'Importance': perm_importance.importances_mean,
    'Importance (%)': (perm_importance.importances_mean / total_importance) * 100
}).sort_values(by='Importance', ascending=False)

feature_importances['Original Feature'] = feature_importances['Feature'].apply(lambda x: x.split('=')[0] if '=' in x else x)


grouped_importances = feature_importances.groupby('Original Feature').sum().reset_index()
grouped_importances = grouped_importances.sort_values(by='Importance (%)', ascending=False)
# Mostrar los resultados
print(feature_importances)

#%%%
# Visualización de la importancia en porcentaje
feature_importances = feature_importances[feature_importances['Importance (%)'] > 0.21]
plt.figure(figsize=(18, 18))
plt.barh(feature_importances['Feature'], feature_importances['Importance (%)'], color='skyblue')
plt.xlabel('Importance (%)', fontsize = 20)
plt.xticks([0,20,40,60],fontsize=15) 
plt.xlim(left=0)
plt.yticks(fontsize = 20)
plt.title('Permutation Feature Importances (in %)')
plt.gca().invert_yaxis()
plt.grid(visible = True, axis = 'x')
plt.show()

# visualizacion de las variables originales
plt.figure(figsize=(18, 18))
plt.barh(grouped_importances['Original Feature'], grouped_importances['Importance (%)'], color='skyblue')
plt.xlabel('Importance (%)', fontsize = 20)
plt.xticks(fontsize=15) 
plt.yticks(fontsize = 20)
plt.gca().invert_yaxis()  # To display the most important features at the top
plt.grid(visible = True, axis = 'x') 
plt.show()

