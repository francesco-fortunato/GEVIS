const sqlite3 = require('sqlite3').verbose();
const fs = require('fs');
const path = require('path');
const { parseString } = require('xml2js');

// Database file path
const dbPath = path.join(__dirname, 'miRNA_database.db');
const xgmmlFilePath = path.join(__dirname, 'mirtarbase_hsa_9.0.xgmml');

// Initialize database connection
const db = new sqlite3.Database(dbPath, (err) => {
  if (err) {
    console.error('Error connecting to database:', err);
  } else {
    console.log('Connected to SQLite database.');
    initializeDatabase();
  }
});

// Function to create the table if it doesn't exist
function initializeDatabase() {
  const createTableQuery = `
    CREATE TABLE IF NOT EXISTS miRNA_targets (
      id INTEGER PRIMARY KEY,
      source_miRNA TEXT NOT NULL,
      target_gene TEXT NOT NULL,
      interaction TEXT,
      datasource TEXT,
      experiments TEXT,
      support_type TEXT
    );
  `;

  db.run(createTableQuery, (err) => {
    if (err) {
      console.error('Error creating table:', err);
    } else {
      console.log('Table created or already exists.');
    }
  });
}

function populateDatabase() {
    const insertQuery = `
      INSERT INTO miRNA_targets (source_miRNA, target_gene, interaction, datasource, experiments, support_type)
      VALUES (?, ?, ?, ?, ?, ?);
    `;
  
    const readStream = fs.createReadStream(xgmmlFilePath, { encoding: 'utf8' });
  
    let buffer = '';
    let nodeMapping = {}; // Mapping of node ID to target gene name
    let recordCount = 0;
  
    // First pass: Read nodes and build a mapping
    readStream.on('data', (chunk) => {
      buffer += chunk;
  
      // Process each node element within the chunk
      let start;
      while ((start = buffer.indexOf('<node ')) !== -1) {
        let end = buffer.indexOf('</node>');
        if (end === -1) break;
  
        end += 7; // Account for the length of `</node>`
        const nodeString = buffer.slice(start, end);
        buffer = buffer.slice(end);
  
        parseString(nodeString, (err, result) => {
          if (err) {
            console.error('Error parsing node:', err);
          } else if (result.node) {
            const node = result.node.$;
            const id = node.id;
            const label = node.label;
  
            // Find the target gene name
            let targetGene = null;
            (result.node.att || []).forEach((att) => {
              if (att.$.name === 'Target Gene') {
                targetGene = att.$.value;
              }
            });
  
            // Map the ID to the target gene name
            if (targetGene) {
              nodeMapping[id] = targetGene;
            }
          }
        });
      }
    });
  
    readStream.on('end', () => {
      console.log("1st step completed");
      // Close the initial read stream
      readStream.close();
  
      // Now we need to read the file again for edges
      const readEdgesStream = fs.createReadStream(xgmmlFilePath, { encoding: 'utf8' });
  
      let edgeBuffer = '';
      readEdgesStream.on('data', (chunk) => {
        edgeBuffer += chunk;
  
        // Process each edge element within the chunk
        let start;
        while ((start = edgeBuffer.indexOf('<edge ')) !== -1) {
          let end = edgeBuffer.indexOf('</edge>');
          if (end === -1) break;
  
          end += 7; // Account for the length of `</edge>`
          const edgeString = edgeBuffer.slice(start, end);
          edgeBuffer = edgeBuffer.slice(end);
  
          parseString(edgeString, (err, result) => {
            if (err) {
              console.error('Error parsing edge:', err);
            } else if (result.edge) {
              const edge = result.edge.$;
              const source = edge.source;
              const targetID = edge.target; // This is the ID we need to map
              const interaction = edge.interaction || null;
  
              const targetGene = nodeMapping[targetID] || null; // Lookup the target gene by ID
              // Extract additional information from "att" tags if present
              let datasource = null, experiments = null, supportType = null;
              (result.edge.att || []).forEach((att) => {
                const name = att.$.name;
                const value = att.$.value;
  
                if (name === 'datasource') datasource = value;
                if (name === 'Experiments') experiments = value;
                if (name === 'Support Type') supportType = value;
              });
  
              // Insert data into the database
              db.run(insertQuery, [source, targetGene, interaction, datasource, experiments, supportType], (err) => {
                if (err) {
                  console.error('Error inserting data:', err);
                } else {
                  recordCount++;
                  if (recordCount % 100 === 0) {
                    console.log(`Inserted ${recordCount} records so far...`);
                  }
                }
              });
            }
          });
        }
      });
  
      readEdgesStream.on('end', () => {
        console.log('Database population complete. Total records inserted:', recordCount);
        db.close();
      });
  
      readEdgesStream.on('error', (err) => {
        console.error('Error reading edges file:', err);
      });
    });
  
    readStream.on('error', (err) => {
      console.error('Error reading file:', err);
    });
  }
  
module.exports = db;
