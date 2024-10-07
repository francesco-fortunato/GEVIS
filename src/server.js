const express = require('express');
const bodyParser = require('body-parser');
const path = require('path');
const session = require('express-session');
const FileStore = require('session-file-store')(session); // Session store for file-based storage
const { v4: uuidv4 } = require('uuid'); // Import uuidv4 from uuid module
const app = express();
const port = 3000;

app.use(session({
  store: new FileStore({ 
    // Optional configuration for the file store
    path: './sessions', // Path to store session files
    ttl: 3600 // Session expiration time in seconds (1 hour)
  }),
  secret: 'your_secret_key', // Replace with your actual secret
  resave: false, // Don't save session if unmodified
  saveUninitialized: false, // Don't save uninitialized sessions
  cookie: { secure: false, maxAge: 3600000 } // Session cookie expires in 1 hour
}));

app.use(express.static('public'));
app.use(bodyParser.json({ limit: '50mb' }));

// Handle POST request to /enrichment
app.post('/enrichment', (req, res) => {
  const data = req.body.data;

  // Generate a unique session ID using uuidv4
  req.session.sessionId = generateSessionId();

  // Store the data in the session
  req.session.data = data;

  // Send the session ID back to the client
  res.json({ sessionId: req.session.sessionId });
});

// Serve the enrichment.html file on a GET request
app.get('/enrichment', (req, res) => {
  res.sendFile(path.join(__dirname, '../public', 'enrichment.html'));
});

////////////////////////////////////////////////////////////////////////////////////////////////////

// Handle POST request to /survival
app.post('/survival', (req, res) => {
  const data = req.body.data;

  // Generate a unique session ID using uuidv4
  req.session.sessionId = generateSessionId();

  // Store the data in the session
  req.session.data = data;

  // Send the session ID back to the client
  res.json({ sessionId: req.session.sessionId });
});

// Serve the survival.html file on a GET request
app.get('/survival', (req, res) => {
  res.sendFile(path.join(__dirname, '../public', 'survival.html'));
});


////////////////////////////////////////////////////////////////////////////////////////////////////


// Handle GET request to retrieve data based on session ID
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

// Function to generate a UUID session ID
function generateSessionId() {
  return uuidv4();
}

app.listen(port, () => {
  console.log(`Server is running at http://localhost:${port}`);
});
