import pandas as pd
import pathlib
import toxindex.utils.chemprop as chemprop
import toxindex.utils.simplecache as simplecache

cachedir = pathlib.Path('cache') / 'function_cache' / 'chemprop_predictions_cache'
cachedir.mkdir(parents=True, exist_ok=True)
pred = simplecache.simple_cache(cachedir)(chemprop.chemprop_predict_all)
inchi = 'InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3'
prediction = pred(inchi)

categories = []
strengths = []
for i in range(len(prediction)):
    prediction_categories = prediction[i]['property']['categories']
    for j in range(len(prediction_categories)):
        categories.append(prediction_categories[j]['category'])
        strengths.append(prediction_categories[j]['strength'])

categories = pd.Series(categories)
strengths = pd.Series(strengths)
numerical_predictions = pd.DataFrame({'category': categories, 'strength': strengths})

means = numerical_predictions.groupby('category')['strength'].mean()
stds = numerical_predictions.groupby('category')['strength'].std()
stats = pd.DataFrame({'mean': means, 'std': stds})
# round to 2 decimal places for LLM
stats = stats.round(2)

llm_input = stats.reset_index().to_dict(orient='records')

# # sanity check on the InChIs
# for pred in prediction:
#     if pred['inchi'] != inchi:
#         print(f"Predicted InChI {pred['inchi']} does not match input InChI {inchi}")
#         break