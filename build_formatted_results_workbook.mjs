import fs from "node:fs/promises";
import path from "node:path";
import { SpreadsheetFile, Workbook } from "@oai/artifact-tool";

const rawArgs = process.argv.slice(2);
const positionalArgs = rawArgs.filter((arg) => !arg.startsWith("--"));
const flagArgs = new Map(
  rawArgs
    .filter((arg) => arg.startsWith("--"))
    .map((arg) => {
      const eqIndex = arg.indexOf("=");
      if (eqIndex === -1) {
        return [arg.slice(2), ""];
      }
      return [arg.slice(2, eqIndex), arg.slice(eqIndex + 1)];
    }),
);

function filterArg(value) {
  return value && value !== "all" ? value : undefined;
}

function listFlag(name) {
  const value = flagArgs.get(name);
  if (!value) {
    return null;
  }
  return new Set(value.split(";").map((item) => item.trim()).filter(Boolean));
}

const inputDir = flagArgs.get("input-dir") ?? "D:/UofT/code/RRRDOR/ITC/result/anchored/plots";
const summaryPath = path.join(inputDir, "anchored_results_summary.csv");
const statusPath = path.join(inputDir, "anchored_complete_scenarios_status.csv");
const outputPath = positionalArgs[0] ?? path.join(inputDir, "anchored_results_formatted.xlsx");
const filters = {
  param: filterArg(positionalArgs[1]),
  event: filterArg(positionalArgs[2]),
  hypothesis: filterArg(positionalArgs[3]),
  n: filterArg(positionalArgs[4]),
};
const scenarioKeyFilter = listFlag("scenario-keys");

const metrics = ["bias", "se", "acc", "coverage", "p"];
const methodOrder = [
  "brm", "CMH.an", "RCMH", "lb.an", "lp.an", "rlp.an",
  "GC", "brm.an", "brm.ad.an", "brm.bc.an", "brm.adbc.an",
];
const groupLabel = "complete_random";

function parseCsv(text) {
  const rows = [];
  let row = [];
  let cell = "";
  let inQuotes = false;

  for (let i = 0; i < text.length; i += 1) {
    const ch = text[i];
    const next = text[i + 1];

    if (ch === '"' && inQuotes && next === '"') {
      cell += '"';
      i += 1;
    } else if (ch === '"') {
      inQuotes = !inQuotes;
    } else if (ch === "," && !inQuotes) {
      row.push(cell);
      cell = "";
    } else if ((ch === "\n" || ch === "\r") && !inQuotes) {
      if (ch === "\r" && next === "\n") {
        i += 1;
      }
      row.push(cell);
      if (row.some((value) => value !== "")) {
        rows.push(row);
      }
      row = [];
      cell = "";
    } else {
      cell += ch;
    }
  }

  if (cell.length > 0 || row.length > 0) {
    row.push(cell);
    rows.push(row);
  }

  const header = rows[0];
  return rows.slice(1).map((values) => Object.fromEntries(
    header.map((name, index) => [name, values[index] ?? ""]),
  ));
}

function num(value) {
  const parsed = Number(value);
  return Number.isFinite(parsed) ? parsed : null;
}

function scenarioKey(row) {
  return [row.param, row.event, row.hypothesis, row.n, row.ess].join("|");
}

function scenarioTitle(row) {
  return {
    left: `${row.param},   n=${row.n},`,
    right: `${row.event}, ${row.hypothesis}, ess=${row.ess}`,
  };
}

function sheetNameSafe(name) {
  return name.replace(/[\\/?*:[\]]/g, "_").slice(0, 31);
}

function setColumnWidths(sheet) {
  const widths = [125, 125, 125, 125, 125, 125, 18, 110, 110, 130, 90, 95, 95];
  widths.forEach((width, index) => {
    sheet.getRangeByIndexes(0, index, 1, 1).format.columnWidthPx = width;
  });
}

function styleTable(sheet, startRow, rowCount) {
  const titleLeft = sheet.getRangeByIndexes(startRow, 0, 1, 1);
  const titleRight = sheet.getRangeByIndexes(startRow, 1, 1, 5);
  const methodHeader = sheet.getRangeByIndexes(startRow + 1, 0, 2, 1);
  const groupHeader = sheet.getRangeByIndexes(startRow + 1, 1, 1, 5);
  const metricHeader = sheet.getRangeByIndexes(startRow + 2, 1, 1, 5);
  const body = sheet.getRangeByIndexes(startRow + 3, 0, rowCount, 6);
  const numericBody = sheet.getRangeByIndexes(startRow + 3, 1, rowCount, 5);

  titleLeft.format = {
    fill: "#FFE699",
    font: { bold: true, color: "#000000", name: "Courier New", size: 11 },
    verticalAlignment: "center",
  };
  titleRight.format = {
    font: { bold: true, color: "#000000", name: "Courier New", size: 11 },
    verticalAlignment: "center",
  };

  methodHeader.merge();
  methodHeader.format = {
    fill: "#1F4E79",
    font: { bold: true, color: "#FFFFFF", name: "Courier New", size: 11 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
  };

  groupHeader.merge();
  groupHeader.format = {
    fill: "#305496",
    font: { bold: true, color: "#FFFFFF", name: "Courier New", size: 11 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
  };

  metricHeader.format = {
    fill: "#DDEBF7",
    font: { bold: true, color: "#000000", name: "Courier New", size: 11 },
    horizontalAlignment: "center",
    verticalAlignment: "center",
  };

  body.format = {
    font: { name: "Courier New", size: 10 },
    borders: { preset: "all", style: "thin", color: "#C9DAF8" },
    verticalAlignment: "center",
  };
  numericBody.format = {
    font: { name: "Courier New", size: 10 },
    horizontalAlignment: "right",
    numberFormat: "0.000",
    borders: { preset: "all", style: "thin", color: "#C9DAF8" },
  };
}

function writeFormattedTables(sheet, rows) {
  setColumnWidths(sheet);
  sheet.showGridLines = false;

  const groups = new Map();
  for (const row of rows) {
    const key = scenarioKey(row);
    if (!groups.has(key)) {
      groups.set(key, []);
    }
    groups.get(key).push(row);
  }

  const sortedGroups = [...groups.values()].sort((a, b) => {
    const av = a[0];
    const bv = b[0];
    return av.param.localeCompare(bv.param)
      || av.event.localeCompare(bv.event)
      || av.hypothesis.localeCompare(bv.hypothesis)
      || Number(av.n) - Number(bv.n)
      || Number(av.ess) - Number(bv.ess);
  });

  let startRow = 0;
  for (const group of sortedGroups) {
    const byMethod = new Map(group.map((row) => [row.method, row]));
    const title = scenarioTitle(group[0]);
    const tableRows = methodOrder
      .filter((method) => byMethod.has(method))
      .map((method) => {
        const row = byMethod.get(method);
        return [method, ...metrics.map((metric) => num(row[metric]))];
      });

    const block = [
      [title.left, title.right, null, null, null, null],
      ["method", groupLabel, null, null, null, null],
      [null, ...metrics],
      ...tableRows,
    ];
    sheet.getRangeByIndexes(startRow, 0, block.length, 6).values = block;
    styleTable(sheet, startRow, tableRows.length);

    tableRows.forEach((values, offset) => {
      const pValue = values[5];
      if (typeof pValue === "number" && pValue < 0.05) {
        sheet.getRangeByIndexes(startRow + 3 + offset, 5, 1, 1).format = {
          fill: "#FFD966",
          font: { name: "Courier New", size: 10 },
          horizontalAlignment: "right",
          numberFormat: "0.000",
          borders: { preset: "all", style: "thin", color: "#C9DAF8" },
        };
      }
    });

    startRow += block.length + 2;
  }

  sheet.freezePanes.freezeRows(0);
}

function writeStatusSheet(sheet, rows) {
  sheet.showGridLines = false;
  sheet.getRange("A1:H1").values = [["param", "event", "hypothesis", "n", "ess", "complete", "missing_R", "missing_RLP"]];
  const values = rows.map((row) => [
    row.param,
    row.event,
    row.hypothesis,
    num(row.n),
    num(row.ess),
    row.complete,
    row.missing_R,
    row.missing_RLP,
  ]);
  sheet.getRangeByIndexes(1, 0, values.length, 8).values = values;
  sheet.getRangeByIndexes(0, 0, 1, 8).format = {
    fill: "#305496",
    font: { bold: true, color: "#FFFFFF" },
  };
  sheet.getRangeByIndexes(0, 0, values.length + 1, 8).format.borders = {
    preset: "all",
    style: "thin",
    color: "#D9EAF7",
  };
  [90, 90, 120, 60, 70, 90, 180, 105].forEach((width, index) => {
    sheet.getRangeByIndexes(0, index, 1, 1).format.columnWidthPx = width;
  });
  sheet.freezePanes.freezeRows(1);
}

const summaryRows = parseCsv(await fs.readFile(summaryPath, "utf8"));
const statusRows = parseCsv(await fs.readFile(statusPath, "utf8"));
const rowMatchesFilters = (row) => {
  if (scenarioKeyFilter && !scenarioKeyFilter.has(scenarioKey(row))) return false;
  if (filters.param && row.param !== filters.param) return false;
  if (filters.event && row.event !== filters.event) return false;
  if (filters.hypothesis && row.hypothesis !== filters.hypothesis) return false;
  if (filters.n && String(row.n) !== String(filters.n)) return false;
  return true;
};
const filteredSummaryRows = summaryRows.filter(rowMatchesFilters);
const filteredStatusRows = statusRows.filter(rowMatchesFilters);

const workbook = Workbook.create();
const formatted = workbook.worksheets.add(sheetNameSafe("Formatted Tables"));
const status = workbook.worksheets.add(sheetNameSafe("Scenario Status"));

writeFormattedTables(formatted, filteredSummaryRows);
writeStatusSheet(status, filteredStatusRows);

await fs.mkdir(inputDir, { recursive: true });

const preview = await workbook.render({
  sheetName: "Formatted Tables",
  range: "A1:F32",
  scale: 1,
  format: "png",
});
await fs.writeFile(
  path.join(inputDir, "anchored_results_formatted_preview.png"),
  new Uint8Array(await preview.arrayBuffer()),
);

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);

console.log(outputPath);
