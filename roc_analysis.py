import pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score
import matplotlib.pyplot as plt

# Path to your data
data_path = cds/.csv

# Load your data
data = pd.read_csv(data_path)

# Assume 'true_labels' are the actual labels and 'predicted_scores' are the scores from your model
true_labels = data['true_labels']
predicted_scores = data['predicted_scores']

# Calculate the ROC curve
fpr, tpr, thresholds = roc_curve(true_labels, predicted_scores)

# Calculate the AUC (Area Under the Curve)
auc = roc_auc_score(true_labels, predicted_scores)
print(f'AUC: {auc}')

# Plot the ROC curve
plt.figure()
plt.plot(fpr, tpr, color='blue', lw=2, label='ROC curve (area = %0.2f)' % auc)
plt.plot([0, 1], [0, 1], color='red', lw=2, linestyle='--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver Operating Characteristic')
plt.legend(loc="lower right")
plt.show()



