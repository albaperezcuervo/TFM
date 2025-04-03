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

# ==============================================
# Required Libraries
# ==============================================

import numpy as np
import pandas as pd
import h5py
import pygad
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, recall_score
from sklearn.ensemble import RandomForestClassifier

from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Input, Dense, Conv2D, MaxPooling2D, Flatten, Dropout, BatchNormalization
from tensorflow.keras.optimizers import Adam
from scikeras.wrappers import KerasClassifier

# ======================================= #  CNN MODEL - Convolutional Neural Network # ==============================================

# Centered Log-Ratio (CLR) normalization to avoid issues with compositional data

def clr_normalization(data):
    data += 1e-6  # Avoid log(0)
    log_data = np.log(data)
    clr_data = log_data - np.mean(log_data, axis=1, keepdims=True)
    return clr_data

# Load microbiota data (tensor format) and metadata labels
with h5py.File('Data/corregido/merged_abundance_acumulative.h5', 'r') as hf:
    rd_array = hf['dataset_name'][:]

metadata = np.load('Data/metadata_array.npy', allow_pickle=True)
labels = metadata[:, 1].astype(int)

# Apply CLR normalization
X_clr = clr_normalization(rd_array)

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(X_clr, labels, test_size=0.2, random_state=42, stratify=labels)

# Flatten data for feature selection
X_train_flat = X_train.reshape(X_train.shape[0], -1)
X_test_flat = X_test.reshape(X_test.shape[0], -1)

# Genetic Algorithm fitness function for CNN

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

# Run Genetic Algorithm to select best features
num_genes = X_train_flat.shape[1]

ga_instance = pygad.GA(
    num_generations=10, num_parents_mating=4,
    fitness_func=fitness_function,
    sol_per_pop=8, num_genes=num_genes,
    gene_space=[0, 1],
    mutation_percent_genes=10
)

ga_instance.run()

# Get best features and reshape for CNN training
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

# Function to build CNN for hyperparameter tuning

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

# Grid search over CNN hyperparameters
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

# Try with cross-validation values of None, 5, 10
cv_results = {}
for cv_val in [None, 5, 10]:
    grid = GridSearchCV(estimator=model, param_grid=param_grid, cv=cv_val, scoring=scoring, refit='accuracy', n_jobs=-1)
    result = grid.fit(X_train_cnn, y_train)
    cv_results[f'cv={cv_val}'] = result

# Train final CNN with best parameters
best_params = cv_results['cv=None'].best_params_
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
    epochs=100,
    validation_data=(X_test_cnn, y_test),
    verbose=1
)

# Evaluate the trained model on test data
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import classification_report, confusion_matrix

# Calculate performance metrics
loss, accuracy, recall = best_cnn.evaluate(X_test_cnn, y_test, verbose=1)
print(f"Test Accuracy: {accuracy:.4f}")
print(f"Test Recall: {recall:.4f}")

# Generate predictions and classification report
predictions = np.round(best_cnn.predict(X_test_cnn))
report = classification_report(y_test, predictions, target_names=['Control (Class 0)', 'Alzheimer (Class 1)'])
print(report)

# Confusion matrix plot
cm = confusion_matrix(y_test, predictions)
plt.figure(figsize=(8, 5))
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=['Control', 'Alzheimer'], yticklabels=['Control', 'Alzheimer'])
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Confusion Matrix - CNN')
plt.tight_layout()
plt.show()

# Smooth curve utility function

def smooth_curve(data, window_size=3):
    smoothed = []
    for i in range(len(data)):
        start = max(0, i - window_size + 1)
        smoothed.append(np.mean(data[start:i + 1]))
    return smoothed

# Plot training vs. validation metrics per epoch (smoothed)
train_accuracy = history.history['accuracy']
val_accuracy = history.history['val_accuracy']
train_recall = history.history['recall']
val_recall = history.history['val_recall']
epochs_range = range(1, len(train_accuracy) + 1)

plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(epochs_range, smooth_curve(train_accuracy), label='Smoothed Train Accuracy')
plt.plot(epochs_range, smooth_curve(val_accuracy), label='Smoothed Validation Accuracy')
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.title('CNN Accuracy over Epochs (Smoothed)')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(epochs_range, smooth_curve(train_recall), label='Smoothed Train Recall', color='green')
plt.plot(epochs_range, smooth_curve(val_recall), label='Smoothed Validation Recall', color='darkgreen')
plt.xlabel('Epochs')
plt.ylabel('Recall')
plt.title('CNN Recall over Epochs (Smoothed)')
plt.legend()

plt.tight_layout()
plt.show()


# ======================================= #  MLPNN MODEL - Multi-Layer Perceptron Network # =============================================

# Load genus-level abundance data and labels
with h5py.File('Data/Merge_Abundance_genus.h5', 'r') as hf:
    rd_array = hf['dataset_name'][:]

metadata = np.load('Data/metadata_array.npy', allow_pickle=True)
labels = metadata[:, 1].astype(int)

# Apply CLR normalization
X_clr = clr_normalization(rd_array)

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(X_clr, labels, test_size=0.2, random_state=42, stratify=labels)

# Genetic Algorithm fitness function for MLP
X_train_flat = X_train
X_test_flat = X_test

def fitness_function(ga_instance, solution, solution_idx):
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

# Run Genetic Algorithm
num_genes = X_train.shape[1]

ga_instance = pygad.GA(
    num_generations=10, num_parents_mating=4,
    fitness_func=fitness_function,
    sol_per_pop=8, num_genes=num_genes,
    gene_space=[0, 1],
    mutation_percent_genes=10
)

ga_instance.run()

# Get best selected features
solution, solution_fitness, _ = ga_instance.best_solution()
selected_features = np.where(solution > 0.5)[0]
X_train_selected = X_train[:, selected_features]
X_test_selected = X_test[:, selected_features]

# Function to build final MLP model for GridSearchCV

def create_mlp(learning_rate=0.001):
    model = Sequential([
        Dense(16, activation='relu', input_shape=(X_train_selected.shape[1],)),
        BatchNormalization(),
        Dense(32, activation='relu'),
        Dropout(0.3),
        Dense(64, activation='relu'),
        BatchNormalization(),
        Dense(128, activation='tanh'),
        Dropout(0.3),
        Dense(128, activation='tanh'),
        BatchNormalization(),
        Dense(256, activation='tanh'),
        Dropout(0.3),
        Dense(1, activation='sigmoid')
    ])
    optimizer = Adam(learning_rate=learning_rate)
    model.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy', 'recall'])
    return model

mlp_model = KerasClassifier(build_fn=create_mlp, epochs=100, batch_size=10, verbose=0)
param_grid = {'model__learning_rate': [0.001, 0.01]}
scoring = {'accuracy': 'accuracy', 'recall': 'recall'}

# GridSearch with different CV strategies
cv_results = {}
for cv_val in [None, 5, 10]:
    grid = GridSearchCV(estimator=mlp_model, param_grid=param_grid, cv=cv_val, scoring=scoring, refit='accuracy', n_jobs=-1)
    result = grid.fit(X_train_selected, y_train)
    cv_results[f'cv={cv_val}'] = result

# Train final MLP with best params
best_params = cv_results['cv=None'].best_params_
best_mlp = create_mlp(learning_rate=best_params['model__learning_rate'])
history = best_mlp.fit(
    X_train_selected, y_train,
    batch_size=10,
    epochs=100,
    validation_data=(X_test_selected, y_test),
    verbose=1
)

# Evaluation and visualization
loss, accuracy, recall = best_mlp.evaluate(X_test_selected, y_test, verbose=1)
print(f"Test Accuracy: {accuracy:.4f}")
print(f"Test Recall: {recall:.4f}")

predictions = np.round(best_mlp.predict(X_test_selected))
report = classification_report(y_test, predictions, target_names=['Control (Class 0)', 'Alzheimer (Class 1)'])
print(report)

cm = confusion_matrix(y_test, predictions)
plt.figure(figsize=(8, 5))
sns.heatmap(cm, annot=True, fmt='d', cmap='Purples', xticklabels=['Control', 'Alzheimer'], yticklabels=['Control', 'Alzheimer'])
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Confusion Matrix - MLP')
plt.tight_layout()
plt.show()

# Smooth training/validation curves
train_accuracy = history.history['accuracy']
val_accuracy = history.history['val_accuracy']
train_recall = history.history['recall']
val_recall = history.history['val_recall']
epochs_range = range(1, len(train_accuracy) + 1)

plt.figure(figsize=(12, 6))
plt.subplot(1, 2, 1)
plt.plot(epochs_range, smooth_curve(train_accuracy), label='Smoothed Train Accuracy')
plt.plot(epochs_range, smooth_curve(val_accuracy), label='Smoothed Validation Accuracy')
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.title('MLP Accuracy over Epochs (Smoothed)')
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(epochs_range, smooth_curve(train_recall), label='Smoothed Train Recall', color='green')
plt.plot(epochs_range, smooth_curve(val_recall), label='Smoothed Validation Recall', color='darkgreen')
plt.xlabel('Epochs')
plt.ylabel('Recall')
plt.title('MLP Recall over Epochs (Smoothed)')
plt.legend()

plt.tight_layout()
plt.show()


# ======================================= #  RANDOM FOREST MODEL # =============================================

# Load abundance data and metadata
with h5py.File('Data/Merge_Abundance_genus.h5', 'r') as hf:
    rd_array = hf['dataset_name'][:]

metadata = np.load('Data/metadata_array.npy', allow_pickle=True)
labels = metadata[:, 1].astype(int)

# Apply CLR normalization
X_clr = clr_normalization(rd_array)

# Train/test split
X_train, X_test, y_train, y_test = train_test_split(X_clr, labels, test_size=0.2, random_state=42, stratify=labels)

X_train_flat = X_train
X_test_flat = X_test

# Genetic Algorithm fitness function for Random Forest

def fitness_function(ga_instance, solution, solution_idx):
    selected_features = np.where(solution > 0.5)[0]
    if len(selected_features) == 0:
        return 0
    X_train_selected = X_train_flat[:, selected_features]
    X_test_selected = X_test_flat[:, selected_features]

    model = RandomForestClassifier(n_estimators=100, max_depth=None, random_state=42)
    model.fit(X_train_selected, y_train)
    y_pred = model.predict(X_test_selected)
    return accuracy_score(y_test, y_pred)

# Run Genetic Algorithm for feature selection
num_genes = X_train.shape[1]

ga_instance = pygad.GA(
    num_generations=10, num_parents_mating=4,
    fitness_func=fitness_function,
    sol_per_pop=8, num_genes=num_genes,
    gene_space=[0, 1],
    mutation_percent_genes=10
)

ga_instance.run()

# Select best features and prepare data
solution, solution_fitness, _ = ga_instance.best_solution()
selected_features = np.where(solution > 0.5)[0]
X_train_selected = X_train[:, selected_features]
X_test_selected = X_test[:, selected_features]

# GridSearchCV for Random Forest hyperparameter tuning
param_grid = {
    'n_estimators': [50, 100, 200],
    'max_depth': [None, 10, 20, 30],
    'min_samples_split': [2, 5, 10],
    'min_samples_leaf': [1, 2, 4]
}

scoring = {'accuracy': 'accuracy', 'recall': 'recall'}
cv_results = {}
for cv_val in [None, 5, 10]:
    grid = GridSearchCV(RandomForestClassifier(random_state=42), param_grid=param_grid, cv=cv_val,
                        scoring=scoring, refit='accuracy', n_jobs=-1)
    result = grid.fit(X_train_selected, y_train)
    cv_results[f'cv={cv_val}'] = result

# Train final Random Forest model
best_params = cv_results['cv=None'].best_params_
best_rf = RandomForestClassifier(
    n_estimators=best_params['n_estimators'],
    max_depth=best_params['max_depth'],
    min_samples_split=best_params['min_samples_split'],
    min_samples_leaf=best_params['min_samples_leaf'],
    random_state=42
)
best_rf.fit(X_train_selected, y_train)

# Evaluate model
y_pred = best_rf.predict(X_test_selected)
test_accuracy = accuracy_score(y_test, y_pred)
test_recall = recall_score(y_test, y_pred)

print(f"Test Accuracy: {test_accuracy:.4f}")
print(f"Test Recall: {test_recall:.4f}")

report = classification_report(y_test, y_pred, target_names=['Control (Class 0)', 'Alzheimer (Class 1)'])
print(report)

# Confusion matrix
cm = confusion_matrix(y_test, y_pred)
plt.figure(figsize=(8, 5))
sns.heatmap(cm, annot=True, fmt='d', cmap='Greens', xticklabels=['Control', 'Alzheimer'], yticklabels=['Control', 'Alzheimer'])
plt.xlabel('Predicted')
plt.ylabel('Actual')
plt.title('Confusion Matrix - Random Forest')
plt.tight_layout()
plt.show()
