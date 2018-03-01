/* rowGetter, label, dataKey, headerRenderer */
import React from 'react';
import {Table, Column, Cell} from 'fixed-data-table-2';
import "fixed-data-table-2/dist/fixed-data-table.css";
import ReactTooltip from 'react-tooltip';
import {findDOMNode} from 'react-dom';


class MyTextCell extends React.Component {
    render() {
        const {rowIndex, field, data} = this.props;
        return (
            <Cell {...this.props}>
                {data[rowIndex][field]}
            </Cell>
       );
    }
}

class MyLinkCell extends React.Component {
    render() {
        const {rowIndex, field, data} = this.props;
        const link = '/' + data[rowIndex][field];
        return (
            <Cell {...this.props}> 
                <a href={link}>{data[rowIndex][field]}</a>
            </Cell>
       );
    }
}

class TooltipCell extends React.Component {
    render() {
        const {data, rowIndex, columnKey} = this.props;
        //const value = data[rowIndex][columnKey];
        const value = columnKey;
        const toolTip = 'Region: ' + data['region'];
        return (
            <Cell
                {...this.props}
                onMouseEnter={() => { ReactTooltip.show(this.refs.valueDiv); }}
                onMouseLeave={() => { ReactTooltip.hide(this.refs.valueDiv); }}>
                <div ref='valueDiv' data-tip={toolTip}>
                {value}
                </div>
            </Cell>
        );
    }
}


class MyTable extends React.Component {

    constructor(props) {
        super(props);

        this.state = {
            rows: [],
            columns: [],
            columnSize: 0,
            filteredDataList: [],
            sortBy: 'ensemble',
            sortDir: 'ASC'};
    }

    componentWillMount() {
        this.jsonList();
    }

    jsonList() {
        fetch('/content/ensemble_list').then(
            function(response){
                return response.json();
            }
        ).then(data => {
                var d = data;
                var columnSize = d.columns.length;

                this.setState({
                    rows: d.data,
                    columns: d.columns,
                    filteredDataList: d.data,
                    columnSize: columnSize
                })
            })
    }

    render() {
        var sortDirArrow = '';
        if (this.state.sortDir !== null){
            sortDirArrow = this.state.sortDir === 'DESC' ? ' ↓' : ' ↑';
        }
        var dataset_columns = []
        for (var index = 0; index < this.state.columnSize; index++) {
            var column_tag = this.state.columns[index]['dataset']
            var column_name = column_tag
            dataset_columns.push(
                <Column 
                    columnKey={column_tag} 
                    //header={<Cell>{column_tag}</Cell>}
                    header={<TooltipCell data={this.state.columns[index]} columnKey={column_tag}></TooltipCell>}
                    cell={<MyTextCell data={this.state.filteredDataList} field={column_tag} />}
                    //cell={<TooltipCell data={this.state.filteredDataList} columnKey={column_tag} />}
                    width={120} 
                    //headerRenderer={this._renderHeader.bind(this,this)}/>)
                />
            )
        }
        return (
            <div>
                <input
                    onChange={this._onFilterChange.bind(this,'ensemble')}
                    placeholder="Filter by Ensemble Name"
                    style={{width:'200px'}}
                />
                <br />
                <Table
                    height={50+((this.state.filteredDataList.length+1) * 30)}
                    width={(this.state.columnSize * 120) + 200}
                    rowsCount={this.state.filteredDataList.length}
                    rowHeight={30}
                    headerHeight={60}>
                    <Column 
                        columnKey="ensemble"
                        header={<Cell>{'Ensemble'+ (this.state.sortBy === 'ensemble' ? sortDirArrow : '')}</Cell>} 
                        cell={<MyLinkCell data={this.state.filteredDataList} field="ensemble" />} 
                        width={200} 
                        fixed={true}
                        //headerRenderer={this._renderHeader.bind(this)}
                    />
                    {dataset_columns}
                </Table>
                <ReactTooltip />
            </div>
        );
    }

    _onFilterChange(cellDataKey, event) {
        if (!event.target.value) {
            this.setState({
                filteredDataList: this.state.rows,
            });
        }
        var filterBy = event.target.value.toString().toLowerCase();
        var size = this.state.rows.length;
        var filteredList = [];
        for (var index = 0; index < size; index++) {
            var v = this.state.rows[index][cellDataKey];
            if (v.toString().toLowerCase().indexOf(filterBy) !== -1) {
                filteredList.push(this.state.rows[index]);
            }
        }
        this.setState({
            filteredDataList: filteredList,
        });
    }

    _sortRowsBy(cellDataKey) {
        var sortDir = this.state.sortDir;
        var sortBy = cellDataKey;
        if (sortBy === this.state.sortBy) {
            sortDir = this.state.sortDir === 'ASC' ? 'DESC' : 'ASC';
        } else {
            sortDir = 'DESC';
        }
        var rows = this.state.filteredDataList.slice();
        rows.sort((a, b) => {
            var sortVal = 0;
            if (a[sortBy] > b[sortBy]) {
                sortVal = 1;
            }
            if (a[sortBy] < b[sortBy]) {
                sortVal = -1;
            }

            if (sortDir === 'DESC') {
                sortVal = sortVal * -1;
            }
            return sortVal;
        });

        this.setState({sortBy, sortDir, filteredDataList : rows});
    }
}

module.exports = MyTable;
