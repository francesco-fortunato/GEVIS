# Use the official Node.js image
FROM node:14

# Create and change to the app directory
WORKDIR /usr/src/app

# Create an empty session folder
RUN mkdir -p sessions

# Copy package.json and package-lock.json
COPY package*.json ./

# Install app dependencies
RUN npm install

# Copy app files
COPY . .

# Expose the port that the app runs on
EXPOSE 3000

# Run the application
CMD ["node", "src/server.js"]
