# -*- coding: utf-8 -*-
"""
Created on Mon Jul  8 09:59:44 2024

@author: santi
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import OneHotEncoder
import numpy as np
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report, roc_curve, roc_auc_score, ConfusionMatrixDisplay, log_loss, brier_score_loss
from sklearn.model_selection import cross_val_predict, GridSearchCV
from sklearn.dummy import DummyClassifier
import joblib

os.chdir('C:/Users/santi/OneDrive/Escritorio/PC-IMBIV/doctorado/asteraceae - muestreo/capitulo 1')


df = pd.read_csv('df_gm.csv')

df.head()


plt.figure(figsize=(20, 15))
sns.histplot(x = "masa", hue = "germina", data = df, bins = 50,  kde = True)
plt.show()

plt.figure(figsize=(20, 15))
sns.boxplot(x = "germina", y = "masa", data = df)
plt.show()

pobs = df['pob'].value_counts().index
dataframes_por_poblacion = {} ## creo diccionario para almacenar las subsdf

for pob in pobs:
    
    df_pob = df[df['pob'] == pob] ## creo subdf
    
    dataframes_por_poblacion[pob] = df_pob



plt.figure(figsize=(40,20))
for i,pob in enumerate(dataframes_por_poblacion) :
    plt.subplot(2, 3, i+1)
    plt.grid(True, linestyle='--', color='gray', alpha=0.5)
    sns.histplot(x = "masa", hue = 'germina', 
                   data = dataframes_por_poblacion[pob], bins = 20,  kde = True)
    plt.title(pob)

plt.show()
    
    
    
## primer modelo de aprendizaje automatico que voy a probar es una regresión binomial
## spliteo







# División entre instancias y etiquetas
X, y = df.loc[:, ['pob', 'masa']], df.germina ## uso solo pob y masa en principio




encoder = OneHotEncoder(handle_unknown='ignore', sparse_output=False)
X_cat = X['pob'].values.reshape(-1,1)

X_cat = encoder.fit_transform(X_cat)

new_columns = []
for col, col_values in zip(['pob'], encoder.categories_):
  for col_value in col_values:
    new_columns.append('{}={}'.format(col, col_value)) ## esto es para darle nombre a la variable codificada por cada valor de la variable original ejemplo Suburb=Albanvale

#  agrego las variables numericas sin modificar
X_masa = X['masa'].values.reshape(-1,1)
X = np.hstack([X_cat, X_masa])
new_columns.append('masa')

X = pd.DataFrame(data=X, columns=new_columns)
print("Matrix has shape {}, with columns: {}".format(X.shape, new_columns))



# división entre entrenamiento y evaluación
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.35, random_state=39)




model = LogisticRegression(random_state=0) ## igual 
model.fit(X_train, y_train) ## en principio la uso por defecto 
y_pred =  model.predict(X_train)
y_pred_test = model.predict(X_test)


print(classification_report(y_train, y_pred))

print(classification_report(y_test, y_pred_test))


confusion_matrix(y_train,y_pred)
confusion_matrix(y_test,y_pred_test)

## funciona bastante bien pero quiero 

## ok voy a usar validacion cruzada

y_pred_cv = cross_val_predict(model, X, y, cv=5)
print(classification_report(y_train, y_pred))

print(classification_report(y_test, y_pred_test))


cm = confusion_matrix(y,y_pred_cv)

disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                              display_labels=model.classes_)
disp.plot()

accuracy_cv = accuracy_score(y, y_pred_cv)
print(f"Precisión de validación cruzada: {accuracy_cv}")
print(classification_report(y,y_pred_cv))
ConfusionMatrixDisplay


## en genral funciona muy bien pero voy a ajustar los hiperparametros 
c = np.logspace(np.log10(1e-05), np.log10(1), num=10).tolist()
prms = [{'penalty':['l2',  'None'], 'C' : c , 'solver': ['lbfgs','liblinear','newton-cholesky', 'saga']},
        {'penalty':['l1'], 'C' : c, 'solver': ['liblinear']},
        {'penalty':['elasticnet'], 'C' : c, 'solver': ['saga']}]

model = LogisticRegression(random_state=0,max_iter = 1000)

model_comp = GridSearchCV(model, prms, scoring = 'neg_log_loss')## uso neg log loss, por que estoy interesado en la presicion de las probabilidades no en las etiquetas.
model_comp.fit(X_train, y_train)


results = model_comp.cv_results_
df_r = pd.DataFrame(results)

df_r = df_r.sort_values(by='rank_test_score', ascending=True)

df_r10 = df_r[df_r['rank_test_score'] <= 20]
df_r10.columns
print('Resultados usando accuracy')
print(df_r10[['param_C', 'param_penalty', 'param_solver', 'mean_test_score', 'std_test_score']])

## me quedo con el modelo numero 3 tiene una penalizacion mas fuerte del modulo de los parametros.

selected_model_params =df_r10['params'][89]

selected_model = LogisticRegression(random_state=0,**selected_model_params)

selected_model.fit(X_train, y_train)

y_pred = selected_model.predict(X_train)
smodel_final_pred = selected_model.predict(X_test)



print(classification_report(y_train, y_pred))

print(classification_report(y_test, smodel_final_pred))


y_pred_prob = selected_model.predict_proba(X_test)[:, 1] # Probabilidad de la clase positiva (1)
log_loss_value = log_loss(y_test, y_pred_prob)
print(f'Log-loss: {log_loss_value}')

# Calcular el Brier Score
brier_score_value = brier_score_loss(y_test, y_pred_prob)
print(f'Brier Score: {brier_score_value}')


# Modelo básico que predice la probabilidad promedio de la clase
dummy_clf = DummyClassifier(strategy='prior') ## va a predecir en base a la distribucion de las clases, calcula la proporción de cada clase en el conjunto de entrenamiento y predice las clases en proporción a esas frecuencias.
dummy_clf.fit(X_train, y_train)

dummy_pred_prob = dummy_clf.predict_proba(X_test)[:,1] # Probabilidad de la clase positiva (1)

# Calcular log-loss y Brier Score del modelo básico
dummy_log_loss = log_loss(y_test, dummy_pred_prob)
dummy_brier_score = brier_score_loss(y_test, dummy_pred_prob)

print(f'Log-loss del modelo básico: {dummy_log_loss}')
print(f'Brier Score del modelo básico: {dummy_brier_score}')
print(f'Log-loss del modelo seleccionado: {log_loss_value}')
print(f'Brier Score del modelo seleccionado: {brier_score_value}')
## veo que el log-loss del modelo es cercano a 0 al igual que el brier score y es muy diferente del modelo dummy, por lo que mi modelo esta funcionando muy bien en el set de testeo. estamos listos para el siguiente paso.
dummy_clf.classes_



## curva roc 

fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
auc_score = roc_auc_score(y_test, y_pred_prob)
print(f"AUC del modelo seleccionado: {auc_score}")

specificity = 1 - fpr
plt.figure()
plt.plot(tpr, specificity, color='red', lw=2, label=f'Curva ROC (área = {auc_score:.2f})')
plt.plot([0, 1], [1, 0], color='blue', linestyle='--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Especificidad')
plt.ylabel('Tasa de Verdaderos Positivos (Sensibilidad)')
plt.title('Curva ROC')
plt.legend(loc="lower right")
plt.show()


cm = confusion_matrix(y_test,smodel_final_pred, labels = [1,0])

disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                              display_labels = np.array([1,0]))
disp.plot()


## ok es un buen modelo
################################################################################################

## voy a probar si el modelo cambia mucho usando solo masa (creo que no)

X, y = df.masa, df.germina ## uso solo pob y masa en principio

X = np.array(X).reshape(-1,1)

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=39)

model = LogisticRegression(random_state=0) ## igual 
model.fit(X_train, y_train) ## en principio la uso por defecto 
y_pred =  model.predict(X_train)
y_pred_test = model.predict(X_test)


print(classification_report(y_train, y_pred))

print(classification_report(y_test, y_pred_test))
## diria que anda exactamente igual 


c = np.logspace(np.log10(1e-05), np.log10(1), num=10).tolist()
prms = [{'penalty':['l2',  'None'], 'C' : c , 'solver': ['lbfgs','liblinear','newton-cholesky', 'saga']},
        {'penalty':['l1'], 'C' : c, 'solver': ['liblinear']},
        {'penalty':['elasticnet'], 'C' : c, 'solver': ['saga']}]

model = LogisticRegression(random_state=0,max_iter = 1000)

model_comp = GridSearchCV(model, prms, scoring = 'neg_log_loss')## uso neg log loss, por que estoy interesado en la presicion de las probabilidades no en las etiquetas.
model_comp.fit(X_train, y_train)


results = model_comp.cv_results_
df_r = pd.DataFrame(results)

df_r = df_r.sort_values(by='rank_test_score', ascending=True)

df_r10 = df_r[df_r['rank_test_score'] <= 20]
df_r10.columns
print('Resultados usando accuracy')
print(df_r10[['param_C', 'param_penalty', 'param_solver', 'mean_test_score', 'std_test_score']])

## modelo [89] {'C': 1.0, 'penalty': 'l1', 'solver': 'liblinear'} como mejor modelo. Dado que el desvio estandar es 0.0594, un coeficiente de varaicion de 0.2 parece ser una variabilidad razonable para la validacion cruzada.Considero que el modelo es lo suficientemente estable y no depende de la particion del set de datos, el comportamiento es predecible independientemente de la paticion.  

## validacion final

selected_model_params =df_r10['params'][89] ## es lo mismo que best_model

selected_model = LogisticRegression(random_state=0,**selected_model_params)

selected_model.fit(X_train, y_train)

y_pred = selected_model.predict(X_train)
smodel_final_pred = selected_model.predict(X_test)



print(classification_report(y_train, y_pred))

print(classification_report(y_test, smodel_final_pred))

## dado que la media de log loss  -0.290256 el desvio estandar es 0.0594, un coeficiente de varaicion de 0.2 parece ser una variabilidad razonable para la validacion cruzada. Considero que el modelo es lo suficientemente estable y no depende de la particion del set de datos, es predecible. 

cm = confusion_matrix(y_test,smodel_final_pred, labels = [1,0])

disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                              display_labels = np.array([1,0]))
disp.plot() ### solo 4 falsos positivos y 2 falsos negativos. 
plt.show()

## comparo con dummy model
y_pred_prob = selected_model.predict_proba(X_test)[:, 1] # Probabilidad de la clase positiva (1)
log_loss_value = log_loss(y_test, y_pred_prob)
print(f'Log-loss: {log_loss_value}')

# Calcular el Brier Score
brier_score_value = brier_score_loss(y_test, y_pred_prob)
print(f'Brier Score: {brier_score_value}')


# Modelo básico que predice la probabilidad promedio de la clase
dummy_clf = DummyClassifier(strategy='prior') ## va a predecir en base a la distribucion de las clases, calcula la proporción de cada clase en el conjunto de entrenamiento y predice las clases en proporción a esas frecuencias.
dummy_clf.fit(X_train, y_train)

dummy_pred_prob = dummy_clf.predict_proba(X_test)[:,1] # Probabilidad de la clase positiva (1)

# Calcular log-loss y Brier Score del modelo básico
dummy_log_loss = log_loss(y_test, dummy_pred_prob)
dummy_brier_score = brier_score_loss(y_test, dummy_pred_prob)

print(f'Log-loss del modelo básico: {dummy_log_loss}')
print(f'Brier Score del modelo básico: {dummy_brier_score}')
print(f'Log-loss del modelo seleccionado: {log_loss_value}')
print(f'Brier Score del modelo seleccionado: {brier_score_value}')


## ok

## curva roc 

fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
auc_score = roc_auc_score(y_test, y_pred_prob)
print(f"AUC del modelo seleccionado: {auc_score}")

specificity = 1 - fpr
plt.figure()
plt.plot(tpr, specificity, color='red', lw=2, label=f'Curva ROC (área = {auc_score:.2f})')
plt.plot([0, 1], [1, 0], color='blue', linestyle='--', lw=2)
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('Especificidad')
plt.ylabel('Tasa de Verdaderos Positivos (Sensibilidad)')
plt.title('Curva ROC')
plt.legend(loc="lower right")
plt.show()
## ultima comprobacion a ver que tan sensibl puede ser
brier_score_values = []
log_loss_values = []
for i in range(1000):
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25)
    s_model = LogisticRegression(**selected_model_params)
    s_model.fit(X_train, y_train)
    y_pred_prob = s_model.predict_proba(X_test)[:, 1]
    brier_score_value = brier_score_loss(y_test, y_pred_prob)
    brier_score_values.append(brier_score_value)
    log_loss_value = log_loss(y_test, y_pred_prob)
    log_loss_values.append(log_loss_value)

plt.figure(figsize=(20,10))
plt.subplot(1,2,1)
plt.plot([x for x in range(1000)],brier_score_values )
plt.yticks(np.arange(0.0, 0.2, step=0.025))
plt.title('brier score')
plt.subplot(1,2,2)
plt.plot([x for x in range(1000)],log_loss_values )
plt.title('log_loss score')
plt.yticks(np.arange(0.0, 0.7, step=0.05))
plt.show()
np.mean(log_loss_values)
np.std(log_loss_values)
np.mean(brier_score_values)
np.std(brier_score_values)


## mi conclusion final es que es un excelente modelo para predecir 
## efecto de la masa usando shap
import shap

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, random_state=39)

best_model = model_comp.best_estimator_
#selected_model = LogisticRegression(random_state=0,**selected_model_params)
best_model.fit(X_train, y_train)
best_model.coef_
best_model.intercept_
## visualizacion con shap y lime
# Crear el objeto explainer

explainer = shap.Explainer(best_model, X_train)

# Calcular los valores SHAP

shap_values = explainer(X_test)

# Visualizar los valores SHAP
shap.summary_plot(shap_values, X_test, show = False)
fig, ax = plt.gcf(), plt.gca()
ax.set_ylabel('Masa', rotation=0)
plt.yticks([]) 
plt.show()
## forma manual de dependence_plot
import seaborn as sns
x = X_test.reshape(74,)
y = shap_values.values.reshape(74,)
d = {'SHAP values': y, 'masa': x}
shap = pd.DataFrame(data = d)
plt.figure(figsize = (20,10))
sns.scatterplot(data = shap, x = 'masa', y = 'SHAP values')




## vemos que hay un impacto positivo de valores altos de masa y uno negativo de valores negativos. Hay un fuerte efecto de la masa.
'''
## LIME
from lime.lime_tabular import LimeTabularExplainer
explainer = LimeTabularExplainer(training_data=np.array(X_train),
                                 feature_names=['masa'],
                                 class_names=['No germina', 'Germina'],
                                 mode='classification')

i = 10
exp = explainer.explain_instance(data_row=X_test[i], 
                                 predict_fn=selected_model.predict_proba)

exp.show_in_notebook(show_all=True)
exp.as_pyplot_figure()
'''## lime tiene sentido para mas variables.
## AJUSTO MODELO FINAL CON TODO EL SET DE DATOS.
best_model = model_comp.best_estimator_
best_model.fit(X, y)

best_model.intercept_
best_model.coef_
print(f'masa en la cual la probabilidad de germinacion es 0.5 {-(best_model.intercept_/best_model.coef_)}')

y_pred_prob = best_model.predict_proba(X)[:, 1]
fpr, tpr, thresholds = roc_curve(y, y_pred_prob)
auc_score = roc_auc_score(y, y_pred_prob) ## da igual 

## el threshold que maximiza especificidad y sencibilidad
# Calcular la métrica de Youden
J = tpr - fpr ## maximizo este valor
optimal_threshold_index = J.argmax()
optimal_threshold = thresholds[optimal_threshold_index]
optimal_threshold ## este es el ubmbral que maximiza la especificidad y la sencibilidad, es el que deberia usarse para establecer la linea de corte y es 0.4028
print('el umbrol optimo es', optimal_threshold)
## guardo el modelo final
joblib.dump(best_model, filename = 'modelo_probgerm.pk1')


