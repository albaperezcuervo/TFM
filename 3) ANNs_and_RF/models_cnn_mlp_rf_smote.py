# ============================================== # SETUP & SEEDING # ==============================================

# Clean session and set random seed between models
# from tensorflow.keras import backend as K
# import gc
# K.clear_session()
# gc.collect()
# import random
# import numpy as np
# import tensorflow as tf
# random.seed(42)
# np.random.seed(42)
# tf.random.set_seed(42)

# ======================================= #  CNN MODEL - Convolutional Neural Network # ==============================================

# Required Libraries
import numpy as np
import pygad
import h5py
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split, GridSearchCV
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Conv2D, MaxPooling2D, Flatten, Dense, Dropout, BatchNormalization, Input
from tensorflow.keras.optimizers import Adam
from scikeras.wrappers import KerasClassifier
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, recall_score
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager

# 1. Load Data
with h5py.File('merged_abundance_acumulative.h5', 'r') as hf:
    rd_array = hf['dataset_name'][:]
metadata = np.load('metadata_array.npy', allow_pickle=True)
labels = metadata[:, 1].astype(int)  # 0 = Control, 1 = Alzheimer

# 2. CLR Normalization
def clr_normalization(data):
    data += 1e-6  # Avoid log(0)
    log_data = np.log(data)
    clr_data = log_data - np.mean(log_data, axis=1, keepdims=True)
    return clr_data
X_clr = clr_normalization(rd_array)

# 3. Apply SMOTE before train-test split
smote = SMOTE(sampling_strategy={0: 750, 1: 750}, random_state=42)
X_resampled, y_resampled = smote.fit_resample(X_clr.reshape(X_clr.shape[0], -1), labels)

# 4. Train-Test Split
X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=42, stratify=y_resampled)

# 5. Flatten data for genetic selection
X_train_flat = X_train
X_test_flat = X_test

# 6. Fitness function for Genetic Algorithm
def fitness_function(ga_instance, solution, solution_idx):
    selected_features = np.where(solution > 0.5)[0]
    if len(selected_features) == 0:
        return 0
    X_train_selected = X_train_flat[:, selected_features]
    X_test_selected = X_test_flat[:, selected_features]
    num_features_selected = len(selected_features)
    num_columns = num_features_selected // 6
    if num_features_selected % 6 != 0:
        num_columns -= 1
    num_features_selected = num_columns * 6
    X_train_selected = X_train_selected[:, :num_features_selected]
    X_test_selected = X_test_selected[:, :num_features_selected]
    X_train_cnn = X_train_selected.reshape(X_train_selected.shape[0], num_columns, 6, 1)
    X_test_cnn = X_test_selected.reshape(X_test_selected.shape[0], num_columns, 6, 1)

    model = Sequential([
        Input(shape=(num_columns, 6, 1)),
        Conv2D(32, (3, 3), activation='relu', padding='same'),
        BatchNormalization(),
        MaxPooling2D((2, 1)),
        Conv2D(64, (3, 3), activation='relu', padding='same'),
        BatchNormalization(),
        MaxPooling2D((2, 1)),
        Flatten(),
        Dense(128, activation='tanh'),
        Dropout(0.5),
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer=Adam(learning_rate=0.0001), loss='binary_crossentropy', metrics=['accuracy'])
    model.fit(X_train_cnn, y_train, epochs=5, batch_size=10, verbose=0, validation_data=(X_test_cnn, y_test))
    _, accuracy = model.evaluate(X_test_cnn, y_test, verbose=0)
    return accuracy

# 7. Genetic Algorithm Configuration and Execution
num_genes = X_train_flat.shape[1]
ga_instance = pygad.GA(
    num_generations=10,
    num_parents_mating=4,
    fitness_func=fitness_function,
    sol_per_pop=8,
    num_genes=num_genes,
    gene_space=[0, 1],
    mutation_percent_genes=10
)
ga_instance.run()

# 8. Final CNN Training with Best Features
solution, solution_fitness, _ = ga_instance.best_solution()
selected_features = np.where(solution > 0.5)[0]
X_train_selected = X_train_flat[:, selected_features]
X_test_selected = X_test_flat[:, selected_features]
num_features_selected = len(selected_features)
num_columns = num_features_selected // 6
if num_features_selected % 6 != 0:
    num_columns -= 1
num_features_selected = num_columns * 6
X_train_selected = X_train_selected[:, :num_features_selected]
X_test_selected = X_test_selected[:, :num_features_selected]
X_train_cnn = X_train_selected.reshape(X_train_selected.shape[0], num_columns, 6, 1)
X_test_cnn = X_test_selected.reshape(X_test_selected.shape[0], num_columns, 6, 1)

# CNN Model Builder for GridSearch
def build_cnn(learning_rate=0.0001, filters_1=32, filters_2=64, kernel_size=(3,3), dropout_rate=0.5):
    model = Sequential([
        Input(shape=(num_columns, 6, 1)),
        Conv2D(filters_1, kernel_size, activation='relu', padding='same'),
        BatchNormalization(),
        MaxPooling2D((2, 1)),
        Conv2D(filters_2, kernel_size, activation='relu', padding='same'),
        BatchNormalization(),
        MaxPooling2D((2, 1)),
        Flatten(),
        Dense(128, activation='tanh'),
        Dropout(dropout_rate),
        Dense(1, activation='sigmoid')
    ])
    optimizer = Adam(learning_rate=learning_rate)
    model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy', 'recall'])
    return model

# GridSearchCV for CNN
model = KerasClassifier(model=build_cnn, verbose=0)
param_grid = {
    'model__learning_rate': [0.0001, 0.001],
    'model__filters_1': [16, 32, 64],
    'model__filters_2': [32, 64, 128],
    'model__kernel_size': [(3,3), (5,5)],
    'model__dropout_rate': [0.3, 0.5],
    'batch_size': [10, 20]
}
scoring = {'accuracy': 'accuracy', 'recall': 'recall'}
grid_1 = GridSearchCV(estimator=model, param_grid=param_grid, cv=None, scoring=scoring, refit='accuracy', n_jobs=-1)
grid_result_1 = grid_1.fit(X_train_cnn, y_train)

# Best CNN Model
best_params = grid_result_1.best_params_
best_cnn = build_cnn(
    learning_rate=best_params['model__learning_rate'],
    filters_1=best_params['model__filters_1'],
    filters_2=best_params['model__filters_2'],
    kernel_size=best_params['model__kernel_size'],
    dropout_rate=best_params['model__dropout_rate']
)
history = best_cnn.fit(
    X_train_cnn, y_train,
    batch_size=best_params['batch_size'],
    epochs=50,
    validation_data=(X_test_cnn, y_test),
    verbose=1
)

# Evaluation
# Evaluation as described in models_cnn_mlp_rf notebooks

# ======================================= #  MLPNN MODEL - Multi-Layer Perceptron Network # =============================================

# Required Libraries for MLP
from tensorflow.keras.wrappers.scikit_learn import KerasClassifier as OldKerasClassifier
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense

# 1. Load Data
with h5py.File('Merge_Abundance_genus.h5', 'r') as hf:
    rd_array = hf['dataset_name'][:]
metadata = np.load('metadata_array.npy', allow_pickle=True)
labels = metadata[:, 1].astype(int)

# 2. CLR Normalization
X_clr = clr_normalization(rd_array)

# 3. Apply SMOTE before train-test split
smote = SMOTE(sampling_strategy={0: 750, 1: 750}, random_state=42)
X_resampled, y_resampled = smote.fit_resample(X_clr, labels)

# 4. Train-Test Split
X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=42, stratify=y_resampled)

# 5. Flatten data for genetic selection
X_train_flat = X_train
X_test_flat = X_test

# 6. Fitness function for Genetic Algorithm

def mlp_fitness_function(ga_instance, solution, solution_idx):
    selected_features = np.where(solution > 0.5)[0]
    if len(selected_features) == 0:
        return 0
    X_train_selected = X_train_flat[:, selected_features]
    X_test_selected = X_test_flat[:, selected_features]

    model = Sequential([
        Dense(8, activation='relu', input_shape=(X_train_selected.shape[1],)),
        Dense(16, activation='relu'),
        Dense(32, activation='tanh'),
        Dense(32, activation='tanh'),
        Dense(64, activation='tanh'),
        Dense(1, activation='sigmoid')
    ])
    model.compile(optimizer=Adam(learning_rate=0.001), loss='binary_crossentropy', metrics=['accuracy'])
    model.fit(X_train_selected, y_train, epochs=5, batch_size=10, verbose=0, validation_data=(X_test_selected, y_test))
    _, accuracy = model.evaluate(X_test_selected, y_test, verbose=0)
    return accuracy

# 7. Run Genetic Algorithm
num_genes = X_train_flat.shape[1]
ga_instance = pygad.GA(
    num_generations=10,
    num_parents_mating=4,
    fitness_func=mlp_fitness_function,
    sol_per_pop=8,
    num_genes=num_genes,
    gene_space=[0, 1],
    mutation_percent_genes=10
)
ga_instance.run()

# 8. Get Selected Features
solution, solution_fitness, _ = ga_instance.best_solution()
selected_features = np.where(solution > 0.5)[0]
X_train_selected = X_train_flat[:, selected_features]
X_test_selected = X_test_flat[:, selected_features]

# 9. Build MLP Model for GridSearch

def create_mlp(learning_rate=0.001):
    model = Sequential([
        Dense(8, activation='relu', input_shape=(X_train_selected.shape[1],)),
        Dense(16, activation='relu'),
        Dense(32, activation='tanh'),
        Dense(32, activation='tanh'),
        Dense(64, activation='tanh'),
        Dense(1, activation='sigmoid')
    ])
    optimizer = Adam(learning_rate=learning_rate)
    model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy', 'recall'])
    return model

mlp_model = OldKerasClassifier(build_fn=create_mlp, epochs=50, batch_size=10, verbose=0)
param_grid = {'model__learning_rate': [0.001, 0.01]}
scoring = {'accuracy': 'accuracy', 'recall': 'recall'}

# 10. GridSearchCV
mlp_grid = GridSearchCV(estimator=mlp_model, param_grid=param_grid, cv=5, scoring=scoring, refit='accuracy', n_jobs=-1)
grid_result = mlp_grid.fit(X_train_selected, y_train)

# 11. Best MLP Model
best_params = grid_result.best_params_
best_mlp = create_mlp(learning_rate=best_params['model__learning_rate'])
history = best_mlp.fit(
    X_train_selected, y_train,
    batch_size=10,
    epochs=50,
    validation_data=(X_test_selected, y_test),
    verbose=1
)

# 12. Evaluation
# Evaluation as described in models_cnn_mlp_rf notebooks

# ======================================= #  RANDOM FOREST MODEL # =============================================

# Required Libraries for Random Forest
from sklearn.ensemble import RandomForestClassifier

# 1. Load Data
with h5py.File('Data/Merge_Abundance_genus.h5', 'r') as hf:
    rd_array = hf['dataset_name'][:]
metadata = np.load('Data/metadata_array.npy', allow_pickle=True)
labels = metadata[:, 1].astype(int)

# 2. CLR Normalization
X_clr = clr_normalization(rd_array)

# 3. Apply SMOTE before train-test split
smote = SMOTE(sampling_strategy={0: 750, 1: 750}, random_state=42)
X_resampled, y_resampled = smote.fit_resample(X_clr, labels)

# 4. Train-Test Split
X_train, X_test, y_train, y_test = train_test_split(X_resampled, y_resampled, test_size=0.2, random_state=42, stratify=y_resampled)

# 5. Flatten data for genetic selection
X_train_flat = X_train
X_test_flat = X_test

# 6. Fitness function for Genetic Algorithm (Optional - if feature selection is required)
# You can integrate pygad here if feature selection is desired.

# 7. Define parameter grid for GridSearch
param_grid = {
    'n_estimators': [50, 100, 200],
    'max_depth': [None, 10, 20, 30],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4]
}

# 8. Define metrics
scoring = {'accuracy': 'accuracy', 'recall': 'recall'}

# 9. Run GridSearchCV
rf_model = RandomForestClassifier(random_state=42)
grid = GridSearchCV(estimator=rf_model, param_grid=param_grid, cv=5, scoring=scoring, refit='accuracy', n_jobs=-1)
grid_result = grid.fit(X_train_flat, y_train)

# 10. Train best model
best_params = grid_result.best_params_
best_rf = RandomForestClassifier(
    n_estimators=best_params['n_estimators'],
    max_depth=best_params['max_depth'],
    min_samples_split=best_params['min_samples_split'],
    min_samples_leaf=best_params['min_samples_leaf'],
    random_state=42
)
best_rf.fit(X_train_flat, y_train)

# 11. Evaluation
# Evaluation as described in models_cnn_mlp_rf notebooks


