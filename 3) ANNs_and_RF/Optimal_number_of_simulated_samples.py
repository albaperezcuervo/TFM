# ============================================== 
# SETUP & SEEDING 
# ==============================================

# Clean session and set random seed between models (optional for reproducibility)
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
#Required Libraries
import numpy as np
import h5py
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_score, recall_score
from imblearn.over_sampling import SMOTE
import pygad
from tensorflow.keras.models import Sequential
from tensorflow.keras.layers import Dense, Dropout, BatchNormalization, Input
from tensorflow.keras.optimizers import Adam
# CLR normalization function
def clr_normalization(data):
    data += 1e-6  # Avoid log(0)
    log_data = np.log(data)
    clr_data = log_data - np.mean(log_data, axis=1, keepdims=True)
    return clr_data

# Load data
with h5py.File('merged_abundance_acumulative.h5', 'r') as hf:
    rd_array = hf['dataset_name'][:]
metadata = np.load('metadata_array.npy', allow_pickle=True)
labels = metadata[:, 1].astype(int)
X_clr = clr_normalization(rd_array)

# Sample sizes to test with SMOTE
sample_sizes = [100, 300, 500, 1000, 1200, 1500, 2000, 2300]
results = []  # Store accuracy, precision, and recall

for num_samples in sample_sizes:
    print(f"Evaluating SMOTE with {num_samples} additional samples per class...")

    # Apply SMOTE
    smote = SMOTE(sampling_strategy={0: num_samples + 145, 1: num_samples + 79}, random_state=42)
    X_resampled, y_resampled = smote.fit_resample(X_clr.reshape(X_clr.shape[0], -1), labels)

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(
        X_resampled, y_resampled, test_size=0.2, random_state=42, stratify=y_resampled
    )

    # Genetic Algorithm for Feature Selection
    num_genes = X_train.shape[1]

    def fitness_function(ga_instance, solution, solution_idx):
        selected_features = np.where(solution > 0.5)[0]
        if len(selected_features) == 0:
            return 0

        X_train_selected = X_train[:, selected_features]
        X_test_selected = X_test[:, selected_features]

        model = Sequential([
            Input(shape=(X_train_selected.shape[1],)),
            Dense(64, activation='relu'),
            BatchNormalization(),
            Dropout(0.5),
            Dense(32, activation='relu'),
            Dense(1, activation='sigmoid')
        ])

        model.compile(optimizer=Adam(learning_rate=0.0001),
                      loss='binary_crossentropy',
                      metrics=['accuracy'])

        model.fit(X_train_selected, y_train, epochs=5, batch_size=20, verbose=0,
                  validation_data=(X_test_selected, y_test))

        _, accuracy = model.evaluate(X_test_selected, y_test, verbose=0)
        return accuracy

    ga_instance = pygad.GA(
        num_generations=10, num_parents_mating=4,
        fitness_func=fitness_function,
        sol_per_pop=8, num_genes=num_genes,
        gene_space=[0, 1],
        mutation_percent_genes=10
    )

    ga_instance.run()
    solution, solution_fitness, _ = ga_instance.best_solution()
    selected_features = np.where(solution > 0.5)[0]

    # CNN with selected features
    X_train_selected = X_train[:, selected_features]
    X_test_selected = X_test[:, selected_features]

    model_cnn = Sequential([
        Input(shape=(X_train_selected.shape[1],)),
        Dense(128, activation='relu'),
        BatchNormalization(),
        Dropout(0.5),
        Dense(64, activation='relu'),
        Dense(1, activation='sigmoid')
    ])

    model_cnn.compile(optimizer=Adam(learning_rate=0.0001),
                      loss='binary_crossentropy',
                      metrics=['accuracy'])

    model_cnn.fit(X_train_selected, y_train, epochs=10, batch_size=20, verbose=0,
                  validation_data=(X_test_selected, y_test))

    # Evaluate model
    _, final_accuracy = model_cnn.evaluate(X_test_selected, y_test, verbose=0)

    # Predict and compute metrics
    y_pred = (model_cnn.predict(X_test_selected) > 0.5).astype(int)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)

    # Save results
    results.append((num_samples, final_accuracy, precision, recall))

    print(f"SMOTE={num_samples} -> Accuracy: {final_accuracy:.4f}, Precision: {precision:.4f}, Recall: {recall:.4f}")

# EVALUATION - CNN 
# Print summary of all results
print("\nSummary of results:")
for num_samples, acc, prec, rec in results:
    print(f"SMOTE={num_samples}: Accuracy={acc:.4f}, Precision={prec:.4f}, Recall={rec:.4f}")

import pandas as pd
import plotly.graph_objects as go

# Load result file
df_results = pd.read_csv("results_smote_cnn.csv")
df_results["Num_SMOTE_Samples"] = df_results["Total_Simulated_Data"]
df_results["Total_Simulated_Data"] = 224 + 2 * df_results["Num_SMOTE_Samples"]

# Create dual-axis plot
fig = go.Figure()

# Add Recall line (left axis, dark green)
fig.add_trace(go.Scatter(x=df_results["Total_Simulated_Data"],
                         y=df_results["Recall"],
                         mode='lines',
                         name="Recall",
                         line=dict(color='darkgreen')))

# Add Accuracy line (right axis, red)
fig.add_trace(go.Scatter(x=df_results["Total_Simulated_Data"],
                         y=df_results["Accuracy"],
                         mode='lines',
                         name="Accuracy",
                         line=dict(color='red'),
                         yaxis="y2"))

# Layout configuration
fig.update_layout(
    title=dict(text="Recall & Accuracy vs Total Simulated Data - CNN", font=dict(size=33, family="Arial", weight="bold")),
    xaxis=dict(title=dict(text="Total Simulated Data", font=dict(size=28, family="Arial", weight="bold")),
               tickfont=dict(size=28, family="Arial", weight="bold")),
    yaxis=dict(title=dict(text="Recall Score", font=dict(size=28, family="Arial", weight="bold")),
               tickfont=dict(size=28, family="Arial", weight="bold"), color="darkgreen"),
    yaxis2=dict(title=dict(text="Accuracy Score", font=dict(size=28, family="Arial", weight="bold")),
                tickfont=dict(size=28, family="Arial", weight="bold"), color="red", overlaying="y", side="right"),
    legend=dict(x=0.7, y=0.5, font=dict(size=20, family="Arial", weight="bold"))
)

# Show plot
fig.show()


# ======================================= #  MLPNN MODEL - Multi-Layer Perceptron Network # =============================================
#Required Libraries
import numpy as np
import h5py
import pygad
from sklearn.model_selection import train_test_split
from sklearn.metrics import precision_score, recall_score
from imblearn.over_sampling import SMOTE
from keras.models import Sequential
from keras.layers import Dense, BatchNormalization, Dropout, Input
from keras.optimizers import Adam
from keras.callbacks import ReduceLROnPlateau
# CLR normalization function
def clr_normalization(data):
    data += 1e-6  # Avoid log(0)
    log_data = np.log(data)
    clr_data = log_data - np.mean(log_data, axis=1, keepdims=True)
    return clr_data

# Load data
with h5py.File('Merge_Abundance_genus.h5', 'r') as hf:
    rd_array = hf['dataset_name'][:]

metadata = np.load('metadata_array.npy', allow_pickle=True)
labels = metadata[:, 1].astype(int)  # 0 = Control, 1 = Alzheimer

# CLR normalization
X_clr = clr_normalization(rd_array)

# Sample sizes to test with SMOTE
sample_sizes = [100, 300, 500, 1000, 1200, 1500, 2000, 2300]
results = []  # Store accuracy, precision, and recall

for num_samples in sample_sizes:
    print(f"Evaluating SMOTE with {num_samples} additional samples per class...")

    # Apply SMOTE
    smote = SMOTE(sampling_strategy={0: num_samples + 145, 1: num_samples + 79}, random_state=42)
    X_resampled, y_resampled = smote.fit_resample(X_clr.reshape(X_clr.shape[0], -1), labels)

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(
        X_resampled, y_resampled, test_size=0.2, random_state=42, stratify=y_resampled
    )

    # Genetic Algorithm for Feature Selection
    num_genes = X_train.shape[1]

    def fitness_function(ga_instance, solution, solution_idx):
        selected_features = np.where(solution > 0.5)[0]
        if len(selected_features) == 0:
            return 0

        X_train_selected = X_train[:, selected_features]
        X_test_selected = X_test[:, selected_features]

        model = Sequential([
            Input(shape=(X_train_selected.shape[1],)),
            Dense(64, activation='relu'),
            BatchNormalization(),
            Dropout(0.5),
            Dense(32, activation='relu'),
            Dense(1, activation='sigmoid')
        ])

        model.compile(optimizer=Adam(learning_rate=0.0001),
                      loss='binary_crossentropy',
                      metrics=['accuracy'])

        model.fit(X_train_selected, y_train, epochs=5, batch_size=20, verbose=0,
                  validation_data=(X_test_selected, y_test))

        _, accuracy = model.evaluate(X_test_selected, y_test, verbose=0)
        return accuracy

    ga_instance = pygad.GA(
        num_generations=10, num_parents_mating=4,
        fitness_func=fitness_function,
        sol_per_pop=8, num_genes=num_genes,
        gene_space=[0, 1],
        mutation_percent_genes=10
    )

    ga_instance.run()
    solution, solution_fitness, _ = ga_instance.best_solution()
    selected_features = np.where(solution > 0.5)[0]

    # Apply MLP with selected features
    X_train_selected = X_train[:, selected_features]
    X_test_selected = X_test[:, selected_features]

    model_mlp = Sequential([
        Dense(8, activation='relu', input_shape=(X_train_selected.shape[1],)),
        Dense(16, activation='relu'),
        Dense(32, activation='tanh'),
        Dense(32, activation='tanh'),
        Dense(64, activation='tanh'),
        Dense(1, activation='sigmoid')
    ])

    optimizer = Adam(learning_rate=0.001)
    model_mlp.compile(optimizer=optimizer, loss='binary_crossentropy', metrics=['accuracy'])

    reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.2, patience=5, min_lr=0.001)

    # Train MLP
    history = model_mlp.fit(
        X_train_selected, y_train,
        validation_split=0.2,
        epochs=100,
        shuffle=True,
        callbacks=[reduce_lr],
        verbose=1
    )

    # Evaluate model
    _, final_accuracy = model_mlp.evaluate(X_test_selected, y_test, verbose=0)

    # Predict and compute metrics
    y_pred = (model_mlp.predict(X_test_selected) > 0.5).astype(int)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)

    # Save results
    results.append((num_samples, final_accuracy, precision, recall))

    print(f"SMOTE={num_samples} -> Accuracy: {final_accuracy:.4f}, Precision: {precision:.4f}, Recall: {recall:.4f}")

# Evaluation follows the same steps as described in the CNN section.
# Simply adapt the filenames if needed (e.g., results_smote_MLP.csv).


# ======================================= #  RANDOM FOREST MODEL # =============================================
# Required librarys
rom sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import accuracy_score, precision_score, recall_score
import numpy as np
import h5py
from imblearn.over_sampling import SMOTE
from sklearn.model_selection import train_test_split
import pygad
# CLR normalization function
def clr_normalization(data):
    data += 1e-6  # Avoid log(0)
    log_data = np.log(data)
    clr_data = log_data - np.mean(log_data, axis=1, keepdims=True)
    return clr_data

# Load data
with h5py.File('Data/Merge_Abundance_genus.h5', 'r') as hf:
    rd_array = hf['dataset_name'][:]

metadata = np.load('metadata_array.npy', allow_pickle=True)
labels = metadata[:, 1].astype(int)  # 0 = Control, 1 = Alzheimer

# CLR normalization
X_clr = clr_normalization(rd_array)

# Sample sizes to test with SMOTE
sample_sizes = [100, 300, 500, 1000, 1200, 1500, 2000, 2300]
results = []  # Store accuracy, precision, and recall

for num_samples in sample_sizes:
    print(f"Evaluating SMOTE with {num_samples} additional samples per class...")

    # Apply SMOTE
    smote = SMOTE(sampling_strategy={0: num_samples + 145, 1: num_samples + 79}, random_state=42)
    X_resampled, y_resampled = smote.fit_resample(X_clr.reshape(X_clr.shape[0], -1), labels)

    # Train-test split
    X_train, X_test, y_train, y_test = train_test_split(
        X_resampled, y_resampled, test_size=0.2, random_state=42, stratify=y_resampled
    )

    # Genetic Algorithm for Feature Selection
    num_genes = X_train.shape[1]

    def fitness_function(ga_instance, solution, solution_idx):
        selected_features = np.where(solution > 0.5)[0]
        if len(selected_features) == 0:
            return 0

        X_train_selected = X_train[:, selected_features]
        X_test_selected = X_test[:, selected_features]

        model = RandomForestClassifier(
            bootstrap=True, max_depth=7, max_features='sqrt',
            min_samples_leaf=3, min_samples_split=15,
            n_estimators=100, class_weight='balanced', random_state=42
        )
        model.fit(X_train_selected, y_train)
        y_pred = model.predict(X_test_selected)

        return accuracy_score(y_test, y_pred)

    ga_instance = pygad.GA(
        num_generations=10, num_parents_mating=4,
        fitness_func=fitness_function,
        sol_per_pop=8, num_genes=num_genes,
        gene_space=[0, 1],
        mutation_percent_genes=10
    )

    ga_instance.run()
    solution, solution_fitness, _ = ga_instance.best_solution()
    selected_features = np.where(solution > 0.5)[0]

    # Train final Random Forest with selected features
    X_train_selected = X_train[:, selected_features]
    X_test_selected = X_test[:, selected_features]

    rf = RandomForestClassifier(
        bootstrap=True, max_depth=7, max_features='sqrt',
        min_samples_leaf=3, min_samples_split=15,
        n_estimators=100, class_weight='balanced', random_state=42
    )
    rf.fit(X_train_selected, y_train)
    y_pred = rf.predict(X_test_selected)

    # Compute metrics
    accuracy = accuracy_score(y_test, y_pred)
    precision = precision_score(y_test, y_pred)
    recall = recall_score(y_test, y_pred)

    results.append((num_samples, accuracy, precision, recall))
    print(f"SMOTE={num_samples} -> Accuracy: {accuracy:.4f}, Precision: {precision:.4f}, Recall: {recall:.4f}")


# Evaluation follows the same steps as described in the CNN section.
# Simply adapt the filenames if needed (e.g., results_smote_RF.csv).

