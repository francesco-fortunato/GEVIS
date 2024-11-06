const express = require('express');
const bodyParser = require('body-parser');
const path = require('path');
const session = require('express-session');
const axios = require('axios'); // Import axios for making HTTP requests
const FileStore = require('session-file-store')(session); // Session store for file-based storage
const { v4: uuidv4 } = require('uuid'); // Import uuidv4 from uuid module
const app = express();
const port = 3000;
const db = require('../db/database.js');

app.use(session({
  store: new FileStore({ 
    path: './sessions', // Path to store session files
    ttl: 3600 // Session expiration time in seconds (1 hour)
  }),
  secret: 'your_secret_key', // Replace with your actual secret
  resave: false, // Don't save session if unmodified
  saveUninitialized: false, // Don't save uninitialized sessions
  cookie: { secure: false, maxAge: 3600000 } // Session cookie expires in 1 hour
}));

app.use(express.static('static'));
app.use(bodyParser.json({ limit: '50mb' }));

// Function to generate a UUID session ID
function generateSessionId() {
  return uuidv4();
}

// Define an endpoint to query miRNA targets
app.post('/queryTargets', (req, res) => {
  const miRNAArray = req.body.miRNA; // example: ['hsa-miR-1308', 'hsa-miR-924']
  
  if (!miRNAArray || !Array.isArray(miRNAArray) || miRNAArray.length === 0) {
    return res.status(400).json({ error: 'Invalid miRNA input' });
  }

  console.log(miRNAArray);

  // Create placeholders for the SQL query
  const placeholders = miRNAArray.map(() => '?').join(', ');
  
  // SQL query to retrieve target genes for the specified miRNAs
  const query = `SELECT * FROM miRNA_targets WHERE source_miRNA IN (${placeholders})`;

  db.all(query, miRNAArray, (err, rows) => {
    if (err) {
      console.error('Error querying database:', err);
      res.status(500).json({ error: 'Database query error' });
    } else {
      res.json({ targets: rows });
    }
  });
});

app.get('/retrieveData', (req, res) => {
  const sessionId = req.query.sessionId;

  // Check if session ID exists in session data
  if (req.session && req.session.sessionId === sessionId) {
    const data = req.session.data;
    if (data) {
      // Respond with the stored data
      res.json(data);
    } else {
      res.status(404).json({ error: 'No data found in session' });
    }
  } else {
    res.status(404).json({ error: 'Session ID not found or expired' });
  }
});

app.post('/enrichment', (req, res) => {
  const data = req.body.data;

  // Generate a unique session ID using uuidv4
  req.session.sessionId = generateSessionId();

  // Store the data in the session
  req.session.data = data;

  // Send the session ID back to the client
  res.json({ sessionId: req.session.sessionId });
});

app.post('/network', async (req, res) => {
  const data = req.body.data; // Access the incoming data

  // Prepare the API call to mirnet
  const apiUrl = "http://api.mirnet.ca/table/mir";
  
  try {
      // Make the API call using axios
      const response = await axios.post(apiUrl, {
          org: data.org,
          idOpt: data.idOpt,
          selSource: data.selSource,
          targetOpt: data.targetOpt,
          myList: data.myList
      }, {
          headers: {
              "Content-Type": "application/json"
          }
      });

      // Send the response from the API back to the client
      res.json({
          sessionId: req.session.id,
          apiResponse: response.data // Include the API response in the JSON sent back
      });
  } catch (error) {
      // Handle errors from the API call
      console.error("Error calling the API:", error);
      res.status(500).json({ error: "Error calling the external API." });
  }
});

app.get('/enrichment', (req, res) => {
  res.sendFile(path.join(__dirname, '../static', 'enrichment.html'));
});

app.post('/survival', (req, res) => {
  const data = req.body.data;

  // Generate a unique session ID using uuidv4
  req.session.sessionId = generateSessionId();

  // Store the data in the session
  req.session.data = data;

  // Send the session ID back to the client
  res.json({ sessionId: req.session.sessionId });
});

app.get('/survival', (req, res) => {
  res.sendFile(path.join(__dirname, '../static', 'survival.html'));
});

app.listen(port, () => {
  console.log(`Server is running at http://localhost:${port}`);
});
