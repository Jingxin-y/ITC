import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const sourceCsv = "D:/UofT/code/RRRDOR/ITC/result/anchored/sample_strata_grid/ITC_sample_strata_anchored_grid_summary.csv";
const outputDir = "D:/codexresults/ITC/rorqual/compareITC/anchored/outputs/sample_strata_tables";
const outputXlsx = path.join(outputDir, "anchored_sample_strata_scenario_tables.xlsx");

function parseCsv(text) {
  const rows = [];
  let row = [];
  let cell = "";
  let inQuotes = false;

  for (let i = 0; i < text.length; i++) {
    const ch = text[i];
    const next = text[i + 1];

    if (inQuotes) {
      if (ch === '"' && next === '"') {
        cell += '"';
        i++;
      } else if (ch === '"') {
        inQuotes = false;
      } else {
        cell += ch;
      }
    } else if (ch === '"') {
      inQuotes = true;
    } else if (ch === ",") {
      row.push(cell);
      cell = "";
    } else if (ch === "\n") {
      row.push(cell);
      rows.push(row);
      row = [];
      cell = "";
    } else if (ch !== "\r") {
      cell += ch;
    }
  }

  if (cell.length > 0 || row.length > 0) {
    row.push(cell);
    rows.push(row);
  }
  return rows.filter((r) => r.some((v) => v !== ""));
}

function toNumber(value) {
  if (value === "" || value == null || value === "NA") return null;
  const out = Number(value);
  return Number.isFinite(out) ? out : value;
}

function scenarioKey(row) {
  return [
    row.scenario_n,
    row.event,
    row.hypothesis,
    row.ess_ratio,
  ].join("|");
}

function stratumLabel(row) {
  return `${row.z1}${row.z2}${row.z3}`;
}

function scenarioTitle(scenario) {
  return `n=${scenario.scenario_n}, event=${scenario.event}, hypothesis=${scenario.hypothesis}, ess=${scenario.ess_ratio}`;
}

function scenarioSort(a, b) {
  const nDiff = Number(a.scenario_n) - Number(b.scenario_n);
  if (nDiff !== 0) return nDiff;
  const eventDiff = String(a.event).localeCompare(String(b.event));
  if (eventDiff !== 0) return eventDiff;
  const hypDiff = String(a.hypothesis).localeCompare(String(b.hypothesis));
  if (hypDiff !== 0) return hypDiff;
  return Number(a.ess_ratio) - Number(b.ess_ratio);
}

function buildRowsForScenario(rows) {
  const strata = [...new Set(rows.map(stratumLabel))].sort();
  const byGroupStratum = new Map();
  for (const row of rows) {
    byGroupStratum.set(`${row.group}|${stratumLabel(row)}`, row);
  }

  const valueFor = (group, stratum, field, divideByRep = false) => {
    const row = byGroupStratum.get(`${group}|${stratum}`);
    if (!row) return "";
    const raw = row[field];
    if (raw === "" || raw == null) return "";
    const value = Number(raw);
    if (!Number.isFinite(value)) return "";
    if (!divideByRep) return value;
    const reps = Number(row.replications);
    return reps > 0 ? value / reps : "";
  };

  const rowDefs = [
    ["IPD N", "IPD", "mean_n", false],
    ["IPD Y=0", "IPD", "mean_y0", false],
    ["IPD Y=1", "IPD", "mean_y1", false],
    ["IPD reweighted N", "IPD_REWEIGHTED", "mean_weighted_n", false],
    ["IPD reweighted Y=0", "IPD_REWEIGHTED", "total_weighted_y0", true],
    ["IPD reweighted Y=1", "IPD_REWEIGHTED", "total_weighted_y1", true],
    ["IPD reweighted ESS", "IPD_REWEIGHTED", "mean_weighted_ess", false],
    ["AD N", "AD", "mean_n", false],
    ["AD Y=0", "AD", "mean_y0", false],
    ["AD Y=1", "AD", "mean_y1", false],
  ];

  return {
    strata,
    matrix: rowDefs.map(([label, group, field, divideByRep]) => [
      label,
      ...strata.map((s) => valueFor(group, s, field, divideByRep)),
    ]),
  };
}

function toA1(rowZero, colZero) {
  let n = colZero + 1;
  let col = "";
  while (n > 0) {
    const rem = (n - 1) % 26;
    col = String.fromCharCode(65 + rem) + col;
    n = Math.floor((n - 1) / 26);
  }
  return `${col}${rowZero + 1}`;
}

function setBlock(sheet, startRow, startCol, values) {
  const range = sheet.getRangeByIndexes(startRow, startCol, values.length, values[0].length);
  range.values = values;
  return range;
}

const csv = await fs.readFile(sourceCsv, "utf8");
const parsed = parseCsv(csv);
const headers = parsed[0];
const data = parsed.slice(1).map((r) => {
  const obj = {};
  headers.forEach((h, i) => {
    obj[h] = toNumber(r[i]);
  });
  return obj;
});

const scenarios = new Map();
for (const row of data) {
  const key = scenarioKey(row);
  if (!scenarios.has(key)) scenarios.set(key, []);
  scenarios.get(key).push(row);
}
const scenarioList = [...scenarios.values()]
  .map((rows) => ({ meta: rows[0], rows }))
  .sort((a, b) => scenarioSort(a.meta, b.meta));

const workbook = Workbook.create();
const tablesSheet = workbook.worksheets.add("Scenario Tables");
const longSheet = workbook.worksheets.add("Long Format");
tablesSheet.showGridLines = false;
longSheet.showGridLines = false;

tablesSheet.getRange("A1:I1").merge();
tablesSheet.getRange("A1").values = [["Sample Strata Tables by Scenario"]];
tablesSheet.getRange("A2:I2").merge();
tablesSheet.getRange("A2").values = [[
  "Values are mean counts per simulation replicate unless the row label says ESS. Columns are covariate strata z1z2z3.",
]];
tablesSheet.getRange("A1:I1").format = {
  fill: "#1F4E79",
  font: { bold: true, color: "#FFFFFF", size: 14 },
};
tablesSheet.getRange("A2:I2").format = {
  fill: "#EAF2F8",
  font: { color: "#1F2937", italic: true },
};

let cursor = 4;
for (const scenario of scenarioList) {
  const { strata, matrix } = buildRowsForScenario(scenario.rows);
  const titleRange = tablesSheet.getRangeByIndexes(cursor, 0, 1, strata.length + 1);
  titleRange.merge();
  titleRange.values = [[scenarioTitle(scenario.meta)]];
  titleRange.format = {
    fill: "#D9EAF7",
    font: { bold: true, color: "#17365D" },
  };

  const header = [["Measure", ...strata]];
  const headerRange = setBlock(tablesSheet, cursor + 1, 0, header);
  headerRange.format = {
    fill: "#4472C4",
    font: { bold: true, color: "#FFFFFF" },
    borders: { preset: "all", style: "thin", color: "#FFFFFF" },
  };

  const bodyRange = setBlock(tablesSheet, cursor + 2, 0, matrix);
  bodyRange.format = {
    borders: { preset: "all", style: "thin", color: "#D9E2F3" },
    numberFormat: "0.00",
  };
  tablesSheet.getRangeByIndexes(cursor + 2, 0, matrix.length, 1).format = {
    fill: "#F3F6FA",
    font: { bold: true, color: "#1F2937" },
    borders: { preset: "all", style: "thin", color: "#D9E2F3" },
  };

  cursor += matrix.length + 4;
}

tablesSheet.getRange("A:A").format.columnWidthPx = 185;
tablesSheet.getRange("B:I").format.columnWidthPx = 82;
tablesSheet.freezePanes.freezeRows(2);

const longRows = [[
  "scenario_n", "event", "hypothesis", "ess_ratio", "group",
  "measure", "stratum", "value",
]];
for (const scenario of scenarioList) {
  const { strata, matrix } = buildRowsForScenario(scenario.rows);
  for (const row of matrix) {
    const measure = row[0];
    for (let i = 0; i < strata.length; i++) {
      longRows.push([
        scenario.meta.scenario_n,
        scenario.meta.event,
        scenario.meta.hypothesis,
        scenario.meta.ess_ratio,
        measure.startsWith("IPD reweighted") ? "IPD_REWEIGHTED" :
          measure.startsWith("IPD") ? "IPD" : "AD",
        measure,
        strata[i],
        row[i + 1],
      ]);
    }
  }
}

setBlock(longSheet, 0, 0, longRows);
longSheet.getRangeByIndexes(0, 0, 1, longRows[0].length).format = {
  fill: "#1F4E79",
  font: { bold: true, color: "#FFFFFF" },
};
longSheet.getRangeByIndexes(1, 7, Math.max(longRows.length - 1, 1), 1).format.numberFormat = "0.00";
longSheet.tables.add(`A1:${toA1(longRows.length - 1, longRows[0].length - 1)}`, true, "LongFormatTable");
longSheet.getRange("A:H").format.columnWidthPx = 120;
longSheet.getRange("F:F").format.columnWidthPx = 175;
longSheet.freezePanes.freezeRows(1);

const inspect = await workbook.inspect({
  kind: "table",
  range: "Scenario Tables!A1:I20",
  include: "values",
  tableMaxRows: 20,
  tableMaxCols: 9,
  maxChars: 4000,
});
console.log(inspect.ndjson);

const errors = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 100 },
  summary: "formula error scan",
  maxChars: 2000,
});
console.log(errors.ndjson);

const preview = await workbook.render({
  sheetName: "Scenario Tables",
  range: "A1:I24",
  scale: 1,
  format: "png",
});
await fs.mkdir(outputDir, { recursive: true });
await fs.writeFile(
  path.join(outputDir, "scenario_tables_preview.png"),
  new Uint8Array(await preview.arrayBuffer()),
);

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputXlsx);
console.log(`Saved ${outputXlsx}`);
