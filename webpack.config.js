var webpack = require('webpack');
var path = require('path');
module.exports = {  
  devtool: 'inline-source-map',
  entry: {
      "tabular_ensemble": "./js/containers/ensembleTabularContainer.js",
      "tabular_data_set": "./js/containers/dataSetTabularContainer.js"
  },
  output: {
    path: path.join(__dirname, 'static'),
      filename: "[name].js"
  },
   module: {
    loaders: [
      {
        test: /\.js$/,
        loader: 'babel-loader',
        exclude: /node_modules/,
        query: {
          presets: ['react', 'es2015']
        }
      },
      { 
        test: /\.css$/,
        loader: 'style!css'
      },
      
       { test: /\.jpe?g$|\.gif$|\.png$|\.svg$|\.woff$|\.ttf$|\.wav$|\.mp3$/,
        loader: require.resolve("file-loader") + "?name=../[path][name].[ext]"
       },
      //{ test: /\.jpe?g$|\.gif$|\.png$|\.svg$|\.woff$|\.ttf$|\.eot$/, loader: "url" }
    ]
  },
  plugins: [
  ]
};
