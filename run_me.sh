#!/bin/bash

# Set up and clean the environment
echo "Running flake8..."
flake8 src/

echo "Running unit tests with pytest..."
pytest tests/

echo "Running the main script..."
python src/main.py

echo "Done!"